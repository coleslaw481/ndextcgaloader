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

logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"

DEFAULT_URL = 'https://raw.githubusercontent.com/iVis-at-Bilkent/pathway-mapper/master/samples'

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
        self._networklistfile = args.networklistfile
        self._datadir = os.path.abspath(args.datadir)

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

    def run(self):
        """
        Runs content loading for NDEx TCGA Content Loader
        :param theargs:
        :return:
        """
        self._parse_config()
        self._parse_load_plan()
        for entry in os.listdir(self._datadir):
            fp = os.path.join(self._datadir, entry)
            if not os.path.isfile(fp):
                continue
            self._get_pandas_dataframe(entry)

        return 0

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

        logger.info('lines: ' + str(lines))
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

        logger.info('Edge data frame:\n' + str(edge_df))
        node_df = pd.DataFrame.from_records(node_rows_tuples, columns=node_fields)
        logger.info('Node data frame:\n' + str(node_df))
        df_with_a = edge_df.join(node_df.set_index('NODE_ID'), on='SOURCE')
        df_with_b = edge_df.join(node_df.set_index('NODE_ID'), on='TARGET')
        logger.info('Node data frame A:\n' + str(df_with_a))
        logger.info('Node data frame B:\n' + str(df_with_b))

        df_with_a_b = df_with_a.join(df_with_b.set_index('--NODE_NAME'), on='--NODE_NAME', lsuffix='_A',
                                     rsuffix='_B')
        logger.info('Node data frame A and B:\n' + str(df_with_a_b))
        #df_with_a_b = df_with_a_b.replace('\n', '', regex=True)
        #df_with_a_b['PARTICIPANT_A'] = df_with_a_b['PARTICIPANT_A'].map(lambda x: x.lstrip('[').rstrip(']'))
        #df_with_a_b['PARTICIPANT_B'] = df_with_a_b['PARTICIPANT_B'].map(lambda x: x.lstrip('[').rstrip(']'))
        sys.exit(1)
        return None
        #return df_with_a_b, node_lines, node_fields


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
