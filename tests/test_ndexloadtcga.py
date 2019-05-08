#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ndextcgaloader` package."""

import os
import tempfile
import shutil

import unittest
from ndexutil.config import NDExUtilConfig
from ndextcgaloader import ndexloadtcga
from ndextcgaloader.ndexloadtcga import NDExNdextcgaloaderLoader

import json
import ndex2


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class TestNdextcgaloader(unittest.TestCase):
    """Tests for `ndextcgaloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

        self._package_dir = ndexloadtcga.get_package_dir()
        self._networklistfile = ndexloadtcga.get_networksfile()
        self._loadplan_path = ndexloadtcga.get_load_plan()


        self._networks_for_testing = ndexloadtcga.get_networksdir()
        self._testing_dir = ndexloadtcga.get_testsdir()
        self._sample_networks_in_tests_dir = os.path.join(self._testing_dir, 'sample_networks')

        self._the_args = {
            'conf': None,
            'datadir': ndexloadtcga.get_networksdir(),
            'dataurl': None,
            'loadplan': self._loadplan_path,
            'logconf':  None,
            'networklistfile': self._networklistfile
        }

        self._the_args = dotdict(self._the_args)

        self.NDExTCGALoader = NDExNdextcgaloaderLoader(self._the_args)


    def tearDown(self):
        """Tear down test fixtures, if any."""

    def _test_parse_arguments(self):
        """Tests parse arguments"""
        res = ndexloadtcga._parse_arguments('hi', [])

        self.assertEqual(res.profile, 'ndextcgaloader')
        self.assertEqual(res.verbose, 0)
        self.assertEqual(res.logconf, None)
        self.assertEqual(res.conf, None)

        someargs = ['-vv','--conf', 'foo', '--logconf', 'hi',
                    '--profile', 'myprofy']
        res = ndexloadtcga._parse_arguments('hi', someargs)

        self.assertEqual(res.profile, 'myprofy')
        self.assertEqual(res.verbose, 2)
        self.assertEqual(res.logconf, 'hi')
        self.assertEqual(res.conf, 'foo')


    def _test_setup_logging(self):
        """ Tests logging setup"""
        try:
            ndexloadtcga._setup_logging(None)
            self.fail('Expected AttributeError')
        except AttributeError:
            pass

        # args.logconf is None
        res = ndexloadtcga._parse_arguments('hi', [])
        ndexloadtcga._setup_logging(res)

        # args.logconf set to a file
        try:
            temp_dir = tempfile.mkdtemp()

            logfile = os.path.join(temp_dir, 'log.conf')
            with open(logfile, 'w') as f:
                f.write("""[loggers]
keys=root

[handlers]
keys=stream_handler

[formatters]
keys=formatter

[logger_root]
level=DEBUG
handlers=stream_handler

[handler_stream_handler]
class=StreamHandler
level=DEBUG
formatter=formatter
args=(sys.stderr,)

[formatter_formatter]
format=%(asctime)s %(name)-12s %(levelname)-8s %(message)s""")

            res = ndexloadtcga._parse_arguments('hi', ['--logconf',
                                                                       logfile])
            ndexloadtcga._setup_logging(res)

        finally:
            shutil.rmtree(temp_dir)


    def validate_network(self, network, sample_network, file_name):

        edges, nodes, node_attributes = network.edges, network.nodes, network.nodeAttributes
        sample_edges, sample_nodes, sample_node_attributes = sample_network.edges, sample_network.nodes, sample_network.nodeAttributes


        self.assertEqual(len(edges), len(sample_edges), 'Edges in sample ' + file_name)
        self.assertEqual(len(nodes), len(sample_nodes), 'Nodes in sample ' + file_name)
        self.assertEqual(len(node_attributes), len(sample_node_attributes), 'Node Attributes in sample ' + file_name)

        for key, value in edges.items():
            sample_value = sample_edges[key]
            self.assertEqual(value, sample_value, 'Edges are different in ' + file_name)

        for key, value in nodes.items():
            sample_value = sample_nodes[key]
            self.assertEqual(value, sample_value, 'Nodes are different in ' + file_name)

        for key, value in node_attributes.items():
            sample_value = sample_node_attributes[key]
            self.assertEqual(value, sample_value, 'Nodes attributes are different in ' + file_name)




    def test_main(self):
        """Tests main function"""

        # try where loading config is successful
        try:
            #temp_dir = tempfile.mkdtemp()
            #confile = os.path.join(temp_dir, 'some.conf')
            #with open(confile, 'w') as f:
            #    f.write("""[hi]
            #    {user} = bob
            #    {pw} = smith
            #    {server} = dev.ndexbio.org""".format(user=NDExUtilConfig.USER,
            #                                         pw=NDExUtilConfig.PASSWORD,
            #                                         server=NDExUtilConfig.SERVER))

            #res = ndexloadtcga.main(['myprog.py', '--profile', 'ndextcgaloader'])
            #self.assertEqual(res, 0)

            with open(self._networklistfile, 'r') as networks:
                list_of_network_files = networks.read().splitlines()
                list_of_network_files_in_cx = [network.replace('.txt', '.cx') for network in list_of_network_files]

            with open(self._loadplan_path, 'r') as f:
                self._loadplan = json.load(f)

            count = 1

            self.NDExTCGALoader.parse_load_plan()
            self.NDExTCGALoader.prepare_report_directory()

            for network_file in list_of_network_files:
                network_file_in_cx = network_file.replace('.txt', '.cx')

                # generate NiceCX from network_file
                df, network_description, id_to_gene_dict = self.NDExTCGALoader.get_pandas_dataframe(network_file)

                network = self.NDExTCGALoader.generate_nice_cx_from_panda_df(df,
                                                                network_file, network_description, id_to_gene_dict)

                path_to_sample_network = os.path.join(self._sample_networks_in_tests_dir, network_file_in_cx)

                network_sample_in_cx = ndex2.create_nice_cx_from_file(path_to_sample_network)

                self.validate_network(network, network_sample_in_cx, network_file)

                print('{}) netwok {} passed'.format(count, network_file.replace('.cx', '')))
                count += 1


        finally:
            print('done')




        #finally:
        #    shutil.rmtree(temp_dir)
