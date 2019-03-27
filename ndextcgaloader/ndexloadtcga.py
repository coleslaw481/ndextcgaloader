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

logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"

DEFAULT_URL = 'https://raw.githubusercontent.com/iVis-at-Bilkent/pathway-mapper/master/samples'

# Simple dictionary mapping values in type field to
# normalized values
NODE_TYPE_MAPPING = {'GENE': 'protein',
                     'FAMILY': 'proteinfamily',
                     'COMPLEX': 'complex'}

# name of CX aspect that contains coordinates for nodes
CARTESIANLAYOUT_ASPECT_NAME = 'cartesianLayout'

POSX_NODE_ATTR = 'POSX'
POSY_NODE_ATTR = 'POSY'


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

        for entry in os.listdir(self._datadir):
            if not entry.endswith('.txt'):
                continue
            fp = os.path.join(self._datadir, entry)
            if not os.path.isfile(fp):
                continue
            self._process_file(fp)

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
                logger.debug('No position attribute for node: ' + str(node))
                continue
            posy = network.get_node_attribute(id, POSY_NODE_ATTR)
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

    def _process_file(self, file_name):
        """PRocesses  a file"""
        df, node_lines, node_fields = self._get_pandas_dataframe(file_name)
        if df is None:
            return
        network = t2n.convert_pandas_to_nice_cx_with_load_plan(df, self._loadplan)

        self._remove_nan_nodes(network)
        self._add_coordinates_aspect_from_pos_attributes(network)
        network.set_name(os.path.basename(file_name).replace('.txt', ''))

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

    def _download_data_files(self):
        """
        Downloads data files to temp directory
        :return:
        """

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
        # read file into lines list skipping the starting
        # lines before --NODE_NAME entry
        with open(path_to_file, 'r') as f:
            for line in f:
                if line.startswith('--NODE_NAME'):
                    break
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

        node_df.rename(index=str, columns={'--NODE_NAME': 'NODE'}, inplace=True)
        # node_df['PARENT_ID'] = node_df['PARENT_ID'].map(id_to_gene_dict, na_action='ignore')

        edge_df['SOURCE'] = edge_df['SOURCE'].map(id_to_gene_dict, na_action='ignore')
        edge_df['TARGET'] = edge_df['TARGET'].map(id_to_gene_dict, na_action='ignore')

        edge_df.rename(index=str,
                       columns={'EDGE_TYPEINTERACTION_PUBMED_ID': 'EDGE_TYPE'},
                       inplace=True)

        df_with_a = edge_df.join(node_df.set_index('NODE'), on='SOURCE', how='outer')

        df_with_a_b = df_with_a.join(node_df.set_index('NODE'), on='TARGET',
                                     rsuffix='_B', how='left')
        df_with_a_b = df_with_a_b.astype(str)

        df_with_a_b['NODE_TYPE'] = df_with_a_b['NODE_TYPE'].map(NODE_TYPE_MAPPING, na_action='ignore')
        df_with_a_b['NODE_TYPE_B'] = df_with_a_b['NODE_TYPE_B'].map(NODE_TYPE_MAPPING, na_action='ignore')
        df_with_a_b['EDGE_TYPE'] = df_with_a_b['EDGE_TYPE'].str.lower()

        # for debugging this writes the data frame generated to a file
        # in same directory input tsv files are located
        # with open(path_to_file + '_with_a_b.csv', 'w') as f:
        #    f.write(df_with_a_b.to_csv(sep='\t'))

        return df_with_a_b, node_lines, node_fields


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
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
