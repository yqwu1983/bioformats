#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestVcfGeno2Bed(unittest.TestCase):
    def setUp(self):
        self.__input = os.path.join('data', 'variants', 'test.vcf.gz')
        self.__individuals = os.path.join('data', 'variants',
                                          'individuals.txt')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_vcfgeno2bed(self):
        sys.argv = ['', 'vcfgeno2bed', self.__input, self.__output]
        bioformats.cli.bioformats()
        # check with the specified list of individuals
        sys.argv = ['', 'vcfgeno2bed', self.__input, self.__output,
                    '-i', self.__individuals]
        bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
