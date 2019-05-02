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
                        required=True)
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
    parser.add_argument('--datadir', help='Directory containing data files in '
                                          '--networklistfile', default='./networks')
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

    parser.add_argument('--tcgaversion', help='Version of NDEx TCGA Networks', required=True)

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

    # HGNC Symbol Identifier Pattern defined at https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000362
    HGNC_REGEX = '^[A-Za-z-0-9_]+(\@)?$'

    def __init__(self, args):
        """

        :param args:
        """
        self._conf_file = args.conf
        self._profile = args.profile

        self._tcga_version = args.tcgaversion

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


        self._reportdir = 'reports'

        self._invalid_protein_names_file_path = \
            os.path.join(os.path.abspath(self._reportdir), 'invalid_protein_names.tsv')

        self._nested_nodes_file_path = \
            os.path.join(os.path.abspath(self._reportdir), 'nested_nodes.tsv')

        self._networks_in_cx_dir = 'networks_in_cx'

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


    def _prepare_report_directory(self):
        # create reports directory if it doesn't exist
        if not os.path.exists(self._reportdir):
            os.makedirs(self._reportdir)

        # remove reports (if any) from previous run
        if os.path.exists(self._invalid_protein_names_file_path):
            os.remove(self._invalid_protein_names_file_path)

        if os.path.exists(self._nested_nodes_file_path):
            os.remove(self._nested_nodes_file_path)


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

        self._prepare_report_directory()


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

        networkType=[]
        networkType.append('pathway')
        network.set_network_attribute("type", json.dumps(networkType))

        network.set_network_attribute("version", self._tcga_version)



    def _report_proteins_with_invalid_names(self, node_df, network_name):

        proteins_with_invalid_names = []

        for index, row in node_df.iterrows():
            if (row['NODE_TYPE'] != 'GENE'):
                continue

            protein_name = row['NODE']

            if not re.match(NDExNdextcgaloaderLoader.HGNC_REGEX, protein_name):
                proteins_with_invalid_names.append(protein_name)

        if proteins_with_invalid_names:
            proteins_with_invalid_names.sort()

            with open(self._invalid_protein_names_file_path, 'a+') as f:
                for protein_name in proteins_with_invalid_names:
                    str_to_write = protein_name + '\t' + network_name + '\n'
                    f.write(str_to_write)
                f.write('\n')


    def _get_node_name_and_type(self, node_df, node_id):
        for index, row in node_df.iterrows():
            if (row['NODE_ID'] == node_id):
                return row['NODE'], row['NODE_TYPE']

        return None, None


    def _report_nested_nodes(self, node_df, network_name):

        nested_nodes= []
        nested_nodes_ids = {}

        for index, row in node_df.iterrows():
            if (row['NODE_TYPE'] == 'GENE'):
                continue

            if (row['PARENT_ID'] == '-1'):
                continue

            # we found a node that is not GENE and that has a parent

            parent_node_name, parent_node_type = self._get_node_name_and_type(node_df, row['PARENT_ID'])

            # get name and type of current node
            nested_node_name = row['NODE']
            nested_node_type = row['NODE_TYPE']

            str_to_write = nested_node_name + '\t' + nested_node_type + '\t'
            str_to_write += parent_node_name + '\t' + parent_node_type + '\t' + network_name + '\n'

            nested_nodes.append(str_to_write)
            nested_nodes_ids[row['NODE_ID']] = row['PARENT_ID']

            #if (nested_node_type == 'FAMILY'):
            #    nested_nodes_ids[row['NODE_ID']] = row['PARENT_ID']

        if nested_nodes:

            if not os.path.exists(self._nested_nodes_file_path):
                header = 'nested_node_name\tnested_node_type\tparent_node_name\tparent_node_type\tnetwork\n'
                with open(self._nested_nodes_file_path, 'a+') as f:
                    f.write(header)

            with open(self._nested_nodes_file_path, 'a+') as f:
                f.write('\n')
                for node in nested_nodes:
                    f.write(node)

        return nested_nodes_ids

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

        # now, replace 'name' and 'represents' in network with names;
        # we only have represents for simple nodes (proteins) whose represetns comply with DExNdextcgaloaderLoader.HGNC_REGEX
        for id, node in network.get_nodes():
            node['n'] = id_to_gene_dict[node['n']]

            nodeId = node['@id']

            if network.nodeAttributes and network.nodeAttributes[nodeId]:
                nodeAttributes = network.nodeAttributes[nodeId]

                node_resolvable = False

                for attr in nodeAttributes:

                    if attr['v'] == 'protein':
                        # only simple nodes, i.e. proteins can be  resolvable

                        if re.match(NDExNdextcgaloaderLoader.HGNC_REGEX, id_to_gene_dict[node['r']]):
                            node['r'] = 'hgnc.symbol:' + id_to_gene_dict[node['r']]
                            node_resolvable = True

                        break

                if not node_resolvable:
                    del node['r']

        self._remove_nan_nodes(network)
        self._add_coordinates_aspect_from_pos_attributes(network)
        network.set_name(os.path.basename(file_name).replace('.txt', ''))

        self._set_network_attributes(network, network_description)

        network_update_key = self._net_summaries.get(network.get_name().upper())

        # apply style to network
        network.apply_style_from_network(self._template)

        # save network in CX
        path_to_networks_in_cx = self._networks_in_cx_dir
        if not os.path.exists(path_to_networks_in_cx):
            os.makedirs(path_to_networks_in_cx)

        full_network_in_cx_path = os.path.join(path_to_networks_in_cx, network.get_name() + '.cx')

        with open(full_network_in_cx_path, 'w') as f:
            json.dump(network.to_cx(), f, indent=4)


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

            if bool(re.match(NDExNdextcgaloaderLoader.HGNC_REGEX, element)):
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
        type_complex_or_proteinfamily.append(NODE_TYPE_MAPPING['COMPARTMENT'])

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
                        row['MEMBER'] = '|'.join(sorted(member_node_attributes_set))

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
                        row['MEMBER_B'] = '|'.join(sorted(member_node_attributes_set))

                    member_node_attributes_set.clear()

        return added_parent_id_column_added, added_parent_id_column_b_added


    def _create_names_for_unnamed_nodes(self, df, id_to_gene_dict,
                                        node_type_param, source_param, member_param, node_id_param):
        '''
        :param df:
        :param id_to_gene_dict:
        :param node_type_param: NODE_TYPE or NODE_TYPE_B
        :param source_param: SOURCE or TARGET
        :param member_param: MEMBER or MEMBER_B
        :param node_id_param: NODE_ID or NODE_ID_B
        :return:
        '''

        node_types = ['proteinfamily', 'compartment', 'complex']

        for index, row in df.iterrows():

            node_type = row[node_type_param]

            if node_type not in node_types:
                continue

            source_or_target = row[source_param]
            if source_or_target and source_or_target.strip() and (source_or_target.lower() != 'undefined'):
                continue

            member = row[member_param]
            if not member:
                continue

            # get a list of proteins from member field
            member_proteins = []
            protein_symbols = member.split('|')
            for protein_symbol in protein_symbols:
                protein_array = protein_symbol.split(':')
                if (len(protein_array) > 1):
                    member_proteins.append(protein_array[1])
                else:
                    member_proteins.append(protein_array[0])

            if not member_proteins:
                continue

            member_proteins.sort()

            if (len(member_proteins) > 4):
                member_proteins = member_proteins[0:4]
                member_proteins_str = " ".join(member_proteins)
                member_proteins_str = member_proteins_str + ' ...'
            else:
                member_proteins_str = " ".join(member_proteins)

            node_name = 'family' if (node_type == 'proteinfamily') else node_type

            row[source_or_target] = node_name + ' [ ' + member_proteins_str + ' ]'

            id_to_gene_dict[row[node_id_param]] = row[source_or_target]

        return


    def _normalize_nodes(self, nodes_df, nested_nodes_map):
        for key, value in nested_nodes_map.items():
            nodes_df = nodes_df[nodes_df.NODE_ID != key]
            nodes_df['PARENT_ID'].replace([key], value, inplace=True)

        return nodes_df

    def  _process_nested_nodes(self, node_df, nested_nodes_map, network_name):
        normalized_df_nodes = node_df.copy()

        while nested_nodes_map:
            normalized_df_nodes = self._normalize_nodes(normalized_df_nodes, nested_nodes_map)
            nested_nodes_map = self._report_nested_nodes(normalized_df_nodes, network_name)

        return normalized_df_nodes

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
        network_name = file_name.replace('.txt', '')

        self._report_proteins_with_invalid_names(node_df, network_name)
        nested_nodes_map = self._report_nested_nodes(node_df, network_name)

        if nested_nodes_map:
            node_df = self._process_nested_nodes(node_df, nested_nodes_map, network_name)

        edge_df.rename(index=str,
                       columns={'EDGE_TYPEINTERACTION_PUBMED_ID': 'EDGE_TYPE'},
                       inplace=True)

        # drop duplicate edges: edges that have the same source, target and type need to be removed
        edge_df.drop_duplicates(subset=['SOURCE', 'TARGET', 'EDGE_TYPE'], keep='first', inplace=True)

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

        add_parent_id_column, add_parent_id_column_b = self._add_member_properties(df_with_a_b)
        df_with_a_b = df_with_a_b.replace(np.nan, '', regex=True)


        # now, we remove all nodes that have no edges, or in other words, remove all edges from dataframe that
        # have no Edge Id
        # however, we want to keep nodes that have no edges but
        #    1) (have NODE_TYPE != protein, and PARENT_ID is -1)
        #
        # These nodes represent unrelated to nothing processes (for example, p53/p21 in
        #  BRCA-2012-Cell-cycle-signaling-pathway) and we need to keep them.
        #
        # So, we iterate over df_with_a_b and add rows that satisfy our condition to the new frame, df_final

        df_final = pd.DataFrame()
        for index, row in df_with_a_b.iterrows():
            if (pd.isnull(row['EDGE_ID']) or (row['EDGE_ID']=='')):
                if ((row['NODE_TYPE'] != 'protein') and (row['PARENT_ID'] == '-1')):

                        df_final = df_final.append(row, ignore_index=True)
            else:
                df_final = df_final.append(row, ignore_index=True)


        self._create_names_for_unnamed_nodes(df_final, id_to_gene_dict,'NODE_TYPE', 'SOURCE', 'MEMBER', 'NODE_ID')

        self._create_names_for_unnamed_nodes(df_final, id_to_gene_dict,'NODE_TYPE_B', 'TARGET', 'MEMBER_B','NODE_ID_B')

        # for debugging this writes the data frame generated to a file
        # in same directory input tsv files are located
        path_to_tsv_file = path_to_file.replace('.txt','.tsv')
        with open(path_to_tsv_file, 'w') as f:
            f.write(df_final.to_csv(sep='\t'))

        return df_final, network_description, id_to_gene_dict


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
