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
        self.__vcf = os.path.join('data', 'snpeff', 'snpeff.vcf')
        self.__vcf_no_snpeff = os.path.join('data', 'variants',
                                            'multiallele_indels.vcf')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_vcfgeno2bed(self):
        for input_file in (self.__vcf, self.__vcf_no_snpeff):
            sys.argv = ['', 'vcfeffect2bed', input_file, self.__output]
            bioformats.cli.bioformats()
            sys.argv = ['', 'vcfeffect2bed', input_file, self.__output,
                        '-i', 'HIGH', 'MODERATE']
            bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
