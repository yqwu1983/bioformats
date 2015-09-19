#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


import os
import logging
import sys
import unittest
import tempfile
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestRmOut2Bed(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join(
            'data', 'repeatmasker', 'repeatmasker.out'
        )
        self.__output_file = tempfile.NamedTemporaryFile().name
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_rmout2bed(self):
        # test the default options
        sys.argv = ['', 'rmout2bed', self.__input_file,
                    self.__output_file]
        bioformats.cli.bioformats()
        # test the --short option
        sys.argv = ['', 'rmout2bed', '-s', self.__input_file,
                    self.__output_file]
        bioformats.cli.bioformats()
        # test various values of the --color option
        for color in ('class', 'identity', 'class_identity'):
            sys.argv = ['', 'rmout2bed', '-c', color,
                        self.__input_file, self.__output_file]
            bioformats.cli.bioformats()
        # test various values of the --name option
        for name in ('id', 'name', 'class', 'family', 'class_family'):
            sys.argv = ['', 'rmout2bed', '-n', name,
                        self.__input_file, self.__output_file]
            bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)
