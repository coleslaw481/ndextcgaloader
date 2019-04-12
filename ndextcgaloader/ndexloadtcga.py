#! /usr/bin/env python

import argparse
import sys
import logging
import json
import os
import pandas as pd
from logging import config
from ndexutil.config import NDExUtilConfig
import ndextcgaloader
import ndexutil.tsv.tsv2nicecx2 as t2n
from ndex2.client import Ndex2
import ndex2
import requests

import re

import numpy as np

logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"

DEFAULT_URL = 'https://raw.githubusercontent.com/iVis-at-Bilkent/pathway-mapper/master/samples'

# Simple dictionary mapping values in type field to
# normalized values
NODE_TYPE_MAPPING = {'GENE': 'protein',
                     'FAMILY': 'proteinfamily',
                     'COMPLEX': 'complex',
                     'PROCESS': 'process', # PROCESS is not in vocabulary, so we keep it as is
                     'COMPARTMENT': 'compartment' # COMPARTMENT is not in vocabulary, so we keep it as is
                     }

# name of CX aspect that contains coordinates for nodes
CARTESIANLAYOUT_ASPECT_NAME = 'cartesianLayout'

POSX_NODE_ATTR = 'POSX'
POSY_NODE_ATTR = 'POSY'

POSX_B_NODE_ATTR = 'POSX_B'
POSY_B_NODE_ATTR = 'POSY_B'



def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)
    parser.add_argument('--profile', help='Profile in configuration '
                                          'file to use to load '
                                          'NDEx credentials which means'
                                          'configuration under [XXX] will be'
                                          'used '
                                          '(default '
                                          'ndextcgaloader)',
                        default='ndextcgaloader')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/'
                             'logging.config.html#logging-config-fileformat '
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')

    parser.add_argument('--conf', help='Configuration file to load '
                                       '(default ~/' +
                                       NDExUtilConfig.CONFIG_FILE)
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module and'
                             'in ' + TSV2NICECXMODULE + '. Messages are '
                             'output at these python logging levels '
                             '-v = ERROR, -vv = WARNING, -vvv = INFO, '
                             '-vvvv = DEBUG, -vvvvv = NOTSET (default no '
                             'logging)')
#    parser.add_argument('--dataurl', help='Base URL to use to download networks'
#                                          'listed in --networklistfile file (default ' +
#                        DEFAULT_URL + ')',
#                        default=DEFAULT_URL,
#                        required=True)
#    parser.add_argument('--style', help='Path to NDEx CX file to use for styling'
#                                        'networks', required=True)
    parser.add_argument('--datadir', help='Directory containing data files in '
                                          '--networklistfile')
    parser.add_argument('--loadplan', help='Load plan json file', required=True)
    parser.add_argument('--networklistfile', help='File containing a list of'
                                                  'file names corresponding'
                                                  'to the networks to download'
                                                  'from URL set in --dataurl',
                        required=True)
    parser.add_argument('--style', help='Path to NDEx CX file to use for styling'
                                        'networks', required=True)
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 ndextcgaloader.__version__))

    return parser.parse_args(args)


def _setup_logging(args):
    """
    Sets up logging based on parsed command line arguments.
    If args.logconf is set use that configuration otherwise look
    at args.verbose and set logging for this module and the one
    in ndexutil specified by TSV2NICECXMODULE constant
    :param args: parsed command line arguments from argparse
    :raises AttributeError: If args is None or args.logconf is None
    :return: None
    """

    if args.logconf is None:
        level = (50 - (10 * args.verbose))
        logging.basicConfig(format=LOG_FORMAT,
                            level=level)
        logging.getLogger(TSV2NICECXMODULE).setLevel(level)
        logger.setLevel(level)
        return

    # logconf was set use that file
    logging.config.fileConfig(args.logconf,
                              disable_existing_loggers=False)


class NDExNdextcgaloaderLoader(object):
    """
    Class to load content
    """
    def __init__(self, args):
        """

        :param args:
        """
        self._conf_file = args.conf
        self._profile = args.profile
        self._args = args
        self._user = None
        self._pass = None
        self._server = None
        self._ndex = None
        self._net_summaries = None
        self._networklistfile = args.networklistfile
        self._datadir = os.path.abspath(args.datadir)
        self._template = None
        self._failed_networks = []


    def _parse_config(self):
            """
            Parses config
            :return:
            """
            ncon = NDExUtilConfig(conf_file=self._conf_file)
            con = ncon.get_config()
            self._user = con.get(self._profile, NDExUtilConfig.USER)
            self._pass = con.get(self._profile, NDExUtilConfig.PASSWORD)
            self._server = con.get(self._profile, NDExUtilConfig.SERVER)

    def _parse_load_plan(self):
        """

        :return:
        """
        with open(self._args.loadplan, 'r') as f:
            self._loadplan = json.load(f)

    def _get_user_agent(self):
        """

        :return:
        """
        return 'tcga/' + self._args.version

    def _create_ndex_connection(self):
        """
        creates connection to ndex
        :return:
        """
        if self._ndex is None:
            self._ndex = Ndex2(host=self._server, username=self._user,
                               password=self._pass, user_agent=self._get_user_agent())

    def _load_network_summaries_for_user(self):
        """
        Gets a dictionary of all networks for user account
        <network name upper cased> => <NDEx UUID>
        :return: dict
        """
        net_summaries = self._ndex.get_network_summaries_for_user(self._user)
        self._net_summaries = {}
        for nk in net_summaries:
            if nk.get('name') is not None:
                self._net_summaries[nk.get('name').upper()] = nk.get('externalId')

    def _load_style_template(self):
        """
        Loads the CX network specified by self._args.style into self._template
        :return:
        """
        self._template = ndex2.create_nice_cx_from_file(os.path.abspath(self._args.style))

    def run(self):
        """
        Runs content loading for NDEx TCGA Content Loader
        :param theargs:
        :return:
        """
        self._parse_config()
        self._parse_load_plan()
        self._create_ndex_connection()
        self._load_network_summaries_for_user()
        self._load_style_template()

        with open(self._networklistfile, 'r') as networks:
            list_of_network_files = networks.read().splitlines()
            list_of_network_files.reverse()

        self._download_data_files(DEFAULT_URL, list_of_network_files, self._datadir)

        for network_file in list_of_network_files:
            self._process_file(network_file)

        return 0

    def _remove_edges_to_node(self, network, edge_id_list):
        """
        Removes edges pointing to node
        :param network:
        :param node_id:
        :return:
        """
        for edge_id in edge_id_list:
            network.remove_edge(edge_id)
            e_attrib = network.get_edge_attributes(edge_id)
            if e_attrib is None:
                continue

            for edge_attrs in e_attrib:
                network.remove_edge_attribute(edge_attrs['@id'])

    def _remove_nan_nodes(self, network):
        """
        Removes nodes named nan and any edges to those nodes
        :param network:
        :return: None
        """
        edge_id_list = []
        node_id_list = []
        for id, node in network.get_nodes():
            if node['n'] == 'nan':
                node_id_list.append(id)
                for edge_id, edge_obj in network.get_edges():
                    if edge_obj['s'] == id or edge_obj['t'] == id:
                        edge_id_list.append(edge_id)

        if len(edge_id_list) > 0:
            self._remove_edges_to_node(network, edge_id_list)
        if len(node_id_list) > 0:
            for id in node_id_list:
                network.remove_node(id)

    def _add_coordinates_aspect_from_pos_attributes(self, network):
        """
        Iterates through all nodes in network looking for
        POSX and POSY node attributes. These values are then
        put into the CARTESIAN_LAYOUT aspect. Finally these
        attributes are removed from the nodes
        :return:
        """
        coordlist = []
        for id, node in network.get_nodes():
            posx = network.get_node_attribute(id, POSX_NODE_ATTR)
            if posx is None:
                posx = network.get_node_attribute(id, POSX_B_NODE_ATTR)
                if posx is None:
                    logger.debug('No position attribute for node: ' + str(node))
                    continue

            posy = network.get_node_attribute(id, POSY_NODE_ATTR)

            if posy is None:
                posy = network.get_node_attribute(id, POSY_B_NODE_ATTR)
                if posy is None:
                    logger.debug('No position y attribute for node: ' + str(node) +
                                 ' which is weird cause we got an X position')
                    continue
            coordlist.append({'node': id, 'x': float(posx['v']), 'y': float(posy['v'])})
            network.remove_node_attribute(id, POSX_NODE_ATTR)
            network.remove_node_attribute(id, POSY_NODE_ATTR)
        if len(coordlist) > 0:
            network.set_opaque_aspect(CARTESIANLAYOUT_ASPECT_NAME,
                                      coordlist)

    def _set_network_attributes(self, network, network_description):

        if network_description:
            network.set_network_attribute("description", network_description)

        network.set_network_attribute("prov:wasGeneratedBy", "ndextcgaloader " + ndextcgaloader.__version__)

        #network.set_network_attribute("reference",
        #  'Istemi Bahceci, Ugur Dogrusoz, Konnor C La, Özgün Babur, Jianjiong Gao, Nikolaus Schultz ' +
        #  '<b>PathwayMapper: a collaborative visual web editor for cancer pathways and genomic data.</b><br>' +
        #  '<i>Bioinformatics</i>, Volume 33, Issue 14, 15 July 2017, Pages 2238–2240, ' +
        #  '<a target="_blank" href="https://doi.org/10.1093/bioinformatics/btx149">doi.org/10.1093/bioinformatics/btx149</a>')

        network.set_network_attribute("reference",
          'Francisco Sanchez-Vega, Marco Mina, Joshua Armenia, Walid K.Chatila, Augustin Luna, Konnor C.La, ' +
          'Sofia Dimitriadoy, David L.Liu, Havish S.Kantheti, Sadegh Saghafinia, Debyani Chakravarty, ' +
          'Foysal Daian, Qingsong Gao, Matthew H.Bailey, Wen-Wei Liang, Steven M.Foltz, Ilya Shmulevich, ' +
          'Li Ding, ... Nikolaus Schultz<br> ' +
          '<b>Oncogenic Signaling Pathways in The Cancer Genome Atlas</b><br>' +
          '<i>Cell</i>, Volume 173, Issue 2, 5 April 2018, Pages 321-337.e10 <br>' +
          '<a target="_blank" href="https://doi.org/10.1016/j.cell.2018.03.035">doi: 10.1016/j.cell.2018.03.035</a>')

        network.set_network_attribute("author",
          'Francisco Sanchez-Vega, Marco Mina, Joshua Armenia, Walid K.Chatila, Augustin Luna, Konnor C.La, ' +
          'Sofia Dimitriadoy, David L.Liu, Havish S.Kantheti, Sadegh Saghafinia, Debyani Chakravarty, ' +
          'Foysal Daian, Qingsong Gao, Matthew H.Bailey, Wen-Wei Liang, Steven M.Foltz, Ilya Shmulevich, ' +
          'Li Ding, Nikolaus Schultz')

        network.set_network_attribute("organism", "Human, 9606, Homo sapiens")

    def _process_file(self, file_name):

        id_to_gene_dict = {}
        """Processes  a file"""
        df, network_description, id_to_gene_dict = self._get_pandas_dataframe(file_name)
        if df is None:
            return

        # replace node names with IDs before transforming Panda dataframe to Nice CX;
        # this is done because as of the moment of writing convert_pandas_to_nice_cx_with_load_plan() cannot
        # handle frames with multiple nodes with the same name; so we use unique IDs instead
        for idx, row in df.iterrows():
            row['SOURCE'] = row['NODE_ID']
            row['TARGET'] = row['NODE_ID_B']

        network = t2n.convert_pandas_to_nice_cx_with_load_plan(df, self._loadplan)

        # now, replace 'name' and 'represents' in network with names
        for id, node in network.get_nodes():
            node['n'] = id_to_gene_dict[node['n']]
            node['r'] = id_to_gene_dict[node['r']]

        self._remove_nan_nodes(network)
        self._add_coordinates_aspect_from_pos_attributes(network)
        network.set_name(os.path.basename(file_name).replace('.txt', ''))

        self._set_network_attributes(network, network_description)

        network_update_key = self._net_summaries.get(network.get_name().upper())

        # apply style to network
        network.apply_style_from_network(self._template)

        if network_update_key is not None:
            return network.update_to(network_update_key, self._server, self._user, self._pass,
                                     user_agent=self._get_user_agent())
        else:
            upload_message = network.upload_to(self._server, self._user,
                                               self._pass,
                                               user_agent=self._get_user_agent())
        return upload_message

    def _handle_error(self, network_name):
        print('unable to get network {}'.format(network_name))
        self._failed_networks.append(network_name)

    def _download_data_files(self, tcga_github_repo_url, list_of_networks, output_directory=os.getcwd()):
        """ Downloads data files to temp directory

        This function takes three arguments: URL of repository, list of networks and working directory.
        It downloades all networks specified in list_of_networks from tcga_github_repo_url and
        saves them in output_directory.

        Args:
            tcga_github_repo_url (required): URL of TCGA networks, it is
                https://raw.githubusercontent.com/iVis-at-Bilkent/pathway-mapper/master/samples/

            list_of_networks (required): List of networks in tcga_github_repo_url that need to be downloaded.
                Example of list_of_networks = [
                        'ACC-2016-WNT-signaling-pathway.txt',
                        'BLCA-2014-Histone-modification-pathway.txt',
                        'BLCA-2014-RTK-RAS-PI(3)K-pathway.tx'
                    ]

            output_directory (optional; defaults to current working dir):
                directory on local machine where files from tcga_github_repo_url will be saved,
                Example: '/Users/joe/tcga'

        Returns:
            none
        """

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        for network in list_of_networks:
            try:
                response = requests.get(os.path.join(tcga_github_repo_url, network))

                if response.status_code // 100 == 2:

                    with open(os.path.join(output_directory, network), "w") as received_file:
                        received_file.write(response.content.decode('utf-8-sig'))
                else:
                    self._handle_error(network)

            except requests.exceptions.RequestException as e:
                self._handle_error(network)

        # print list of networks that we failed to download (if any)
        if (self._failed_networks):
            print('failed to receive {} networks:'.format(len(self._failed_networks)))
            for network_name in self._failed_networks:
                print(network_name)

    def _generate_member_node_attributes(self, df):
        l = df.tolist()
        set_of_gene_names_with_prefix = set()
        for element in l:
            if element == 'nan':
                continue

            if bool(re.match('^[A-Za-z-0-9_]+(\@)?$', element)):
                set_of_gene_names_with_prefix.add('hgnc.symbol:' + element)
            else:
                set_of_gene_names_with_prefix.add(element)

        # lisf_of_gene_names_with_prefix = ['hgnc:' + element for element in my_df.tolist()]
        return set_of_gene_names_with_prefix


    def _add_member_properties(self, df):
        added_parent_id_column_added = False
        added_parent_id_column_b_added = False

        type_complex_or_proteinfamily = []
        type_complex_or_proteinfamily.append(NODE_TYPE_MAPPING['FAMILY'])
        type_complex_or_proteinfamily.append(NODE_TYPE_MAPPING['COMPLEX'])

        for idx, row in df.iterrows():
            if row['NODE_TYPE'] in type_complex_or_proteinfamily:
                if not added_parent_id_column_added:
                    added_parent_id_column_added = True
                    df['MEMBER'] = ''
                    break

        if added_parent_id_column_added:
            member_node_attributes_set = set()

            for idx, row in df.iterrows():
                if row['NODE_TYPE'] in type_complex_or_proteinfamily:
                    member_node_id = row['NODE_ID']

                    my_df = df[(df['PARENT_ID'] == member_node_id)]['SOURCE']
                    member_node_attributes = self._generate_member_node_attributes(my_df)
                    member_node_attributes_set.update(member_node_attributes)

                    my_df_b = df[(df['PARENT_ID_B'] == member_node_id)]['TARGET']
                    member_node_attributes = self._generate_member_node_attributes(my_df_b)
                    member_node_attributes_set.update(member_node_attributes)

                    if member_node_attributes_set:
                        mem = '|'.join(member_node_attributes_set)
                        row['MEMBER'] = mem

                    member_node_attributes_set.clear()

        for idx, row in df.iterrows():
            if row['NODE_TYPE_B'] in type_complex_or_proteinfamily:
                if not added_parent_id_column_b_added:
                    added_parent_id_column_b_added = True
                    df['MEMBER_B'] = ''
                    break

        if added_parent_id_column_b_added:
            member_node_attributes_set = set()

            for idx, row in df.iterrows():
                if row['NODE_TYPE_B'] in type_complex_or_proteinfamily:
                    member_node_id = row['NODE_ID_B']

                    # get a list of all target node names with the same id as PARENT_ID_B
                    my_df_b = df[(df['PARENT_ID_B'] == member_node_id)]['TARGET']
                    member_node_attributes = self._generate_member_node_attributes(my_df_b)
                    member_node_attributes_set.update(member_node_attributes)

                    my_df_b = df[(df['PARENT_ID'] == member_node_id)]['SOURCE']
                    member_node_attributes = self._generate_member_node_attributes(my_df_b)
                    member_node_attributes_set.update(member_node_attributes)

                    if member_node_attributes_set:
                        mem = '|'.join(member_node_attributes_set)
                        row['MEMBER_B'] = mem

                    member_node_attributes_set.clear()

        return added_parent_id_column_added, added_parent_id_column_b_added


    def _get_pandas_dataframe(self, file_name):
        """
        Gets pandas data frame from file
        :param file_name:
        :return: tuple (dataframe, node lines list, node fields list)
        """
        path_to_file = os.path.join(os.path.abspath(self._datadir),
                                   file_name)
        if os.path.getsize(path_to_file) is 0:
            logger.error('File is empty: ' + path_to_file)
            return None, None, None

        reached_header = False
        lines = []
        logger.info('Examining file: ' + path_to_file)
        current_line_no = 0;
        network_description = ''
        # read file into lines list skipping the starting
        # lines before --NODE_NAME entry
        with open(path_to_file, 'r') as f:
            for line in f:
                current_line_no += 1
                if line.startswith('--NODE_NAME'):
                    break
                else:
                    if (current_line_no > 1) and len(line.strip()) > 0:
                        network_description += line

            lines.append(line)
            lines.extend(f.readlines())

        mode = "node"
        edge_lines = []
        edge_rows_tuples = []
        node_rows_tuples = []
        node_lines = []
        edge_fields = []
        node_fields = []
        for index in range(len(lines)):
            line = lines[index]
            if index is 0:
                if line.startswith('--'):
                    line = line[2:]
                line = line.strip('\n')
                if line.endswith('--'):
                    line = line[0:-2]
                node_fields = [h.strip() for h in line.split('\t')]
            elif line == '\n':
                mode = "edge_header"
            elif mode is "edge_header":
                edge_fields = [h.strip() for h in line.split('\t')]
                mode = "edge"
            elif mode is "node":
                node_tuple = tuple(line.rstrip().split('\t'))
                node_rows_tuples.append(node_tuple)
                node_lines.append(line)
            elif mode is "edge":
                edge_tuple = tuple(line.rstrip().split('\t'))
                edge_rows_tuples.append(edge_tuple)
                edge_lines.append(line)

        edge_df = pd.DataFrame.from_records(edge_rows_tuples, columns=edge_fields)

        node_df = pd.DataFrame.from_records(node_rows_tuples, columns=node_fields)

        id_to_gene_dict = {}

        for idx, row in node_df.iterrows():
            id_to_gene_dict[row.NODE_ID] = row[0]

        node_df.rename(index=str, columns={'NODE_NAME': 'NODE'}, inplace=True)
        edge_df.rename(index=str, columns={'--EDGE_ID': 'EDGE_ID'}, inplace=True)
        # node_df['PARENT_ID'] = node_df['PARENT_ID'].map(id_to_gene_dict, na_action='ignore')

        # drop duplicate edges: edges that have the same source, target and type need to be removed
        edge_df.drop_duplicates(subset=['SOURCE', 'TARGET', 'EDGE_TYPE'], keep='first', inplace=True)


        # for all nodes where column NODE is empty and NODE_TYPE column is 'FAMILY', set
        # value of NODE to "unnamed family"
        node_df.loc[(node_df['NODE'] == '') & (node_df['NODE_TYPE'] == 'FAMILY'), "NODE"] = "unnamed family"

        edge_df.rename(index=str,
                       columns={'EDGE_TYPEINTERACTION_PUBMED_ID': 'EDGE_TYPE'},
                       inplace=True)

        #edge_df.set_index('EDGE_ID')
        df_with_a = edge_df.join(node_df.set_index('NODE_ID'), on='SOURCE', how='inner')
        df_with_b = edge_df.join(node_df.set_index('NODE_ID'), on='TARGET', how='inner')


        df_with_b = df_with_b.drop(columns=['SOURCE', 'TARGET', 'EDGE_TYPE'])
        df_with_b.rename(index=str, columns={'NODE':'NODE_ID_B'}, inplace=True)
        df_with_b.rename(index=str, columns={'NODE_TYPE':'NODE_TYPE_B'}, inplace=True)
        df_with_b.rename(index=str, columns={'PARENT_ID': 'PARENT_ID_B'}, inplace=True)
        df_with_b.rename(index=str, columns={'POSX': 'POSX_B'}, inplace=True)
        df_with_b.rename(index=str, columns={'POSY': 'POSY_B'}, inplace=True)

        df_with_a_b = df_with_a.join(df_with_b.set_index('EDGE_ID'), on='EDGE_ID', how='right')
        df_with_a_b = df_with_a_b.astype(str)

        df_with_a_b['NODE_TYPE'] = df_with_a_b['NODE_TYPE'].map(NODE_TYPE_MAPPING, na_action='ignore')
        df_with_a_b['NODE_TYPE_B'] = df_with_a_b['NODE_TYPE_B'].map(NODE_TYPE_MAPPING, na_action='ignore')
        df_with_a_b['EDGE_TYPE'] = df_with_a_b['EDGE_TYPE'].str.lower()

        df_with_a_b['NODE'] = df_with_a_b['SOURCE']
        df_with_a_b['NODE_ID_B'] = df_with_a_b['TARGET']

        df_with_a_b['SOURCE'] = edge_df['SOURCE'].map(id_to_gene_dict, na_action='ignore')
        df_with_a_b['TARGET'] = edge_df['TARGET'].map(id_to_gene_dict, na_action='ignore')

        df_with_a_b.rename(index=str, columns={'NODE': 'NODE_ID'}, inplace=True)

        df_with_a_b['NODE_TYPE'].fillna('other', inplace=True)
        #df_with_a_b['NODE_TYPE_B'].fillna('other', inplace=True)


        nodes_with_edges_ids = set()
        for index, row in df_with_a_b.iterrows():
            nodes_with_edges_ids.add(row['NODE_ID'])
            nodes_with_edges_ids.add(row['NODE_ID_B'])

        node_df_without_edges = node_df[~node_df['NODE_ID'].isin(nodes_with_edges_ids)]
        node_df_without_edges['NODE_TYPE'] = node_df_without_edges['NODE_TYPE'].map(NODE_TYPE_MAPPING, na_action='ignore')

        node_df_without_edges.rename(index=str, columns={'NODE': 'SOURCE'}, inplace=True)

        for index, row in node_df_without_edges.iterrows():
            df_with_a_b = df_with_a_b.append(row, ignore_index=True)  # Moving

        df_with_a_b = df_with_a_b.replace(np.nan, '', regex=True)

        add_parent_id_column, add_parent_id_column_b = self._add_member_properties(df_with_a_b)

        # for debugging this writes the data frame generated to a file
        # in same directory input tsv files are located
        with open(path_to_file + '_with_a_b.tsv', 'w') as f:
            f.write(df_with_a_b.to_csv(sep='\t'))

        return df_with_a_b, network_description, id_to_gene_dict


def main(args):
    """
    Main entry point for program
    :param args:
    :return:
    """
    desc = """
    Version {version}

    Loads NDEx TCGA Content Loader data into NDEx (http://ndexbio.org).

    To connect to NDEx server a configuration file must be passed
    into --conf parameter. If --conf is unset the configuration
    the path ~/{confname} is examined.

    The configuration file should be formatted as follows:

    [<value in --profile (default ncipid)>]

    {user} = <NDEx username>
    {password} = <NDEx password>
    {server} = <NDEx server(omit http) ie public.ndexbio.org>


    """.format(confname=NDExUtilConfig.CONFIG_FILE,
               user=NDExUtilConfig.USER,
               password=NDExUtilConfig.PASSWORD,
               server=NDExUtilConfig.SERVER,
               version=ndextcgaloader.__version__)
    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = ndextcgaloader.__version__

    try:
        _setup_logging(theargs)
        loader = NDExNdextcgaloaderLoader(theargs)
        return loader.run()
    except Exception as e:
        logger.exception('Caught exception')
        print('\n\n\tException: {}\n'.format(e))
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
