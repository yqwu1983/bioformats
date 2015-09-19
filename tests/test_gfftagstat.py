#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


import os
import logging
import sys
import unittest
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestGffTagStat(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join('data', 'gff3', 'correct.gff')
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_gfftagstat(self):
        sys.argv = ['', 'gfftagstat', self.__input_file]
        bioformats.cli.bioformats()

        # check filtered input
        sys.argv = ['', 'gfftagstat', self.__input_file, '-s',
                    'example']
        bioformats.cli.bioformats()
