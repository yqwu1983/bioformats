#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.bed
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestVcfGeno2Bed(unittest.TestCase):
    def setUp(self):
        self.__vcf = os.path.join('data', 'variants', 'test.vcf')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_vcf2bed(self):
        sys.argv = ['', 'vcf2bed', self.__vcf, self.__output]
        bioformats.cli.bioformats()
        # make sure that the produced file is a correct BED file
        with open(self.__output) as output_file:
            bed_reader = bioformats.bed.Reader(output_file)
            for record in bed_reader.records():
                self.assertIsInstance(record, bioformats.bed.Record)

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
