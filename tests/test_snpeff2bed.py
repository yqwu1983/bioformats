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


class TestSnpeff2Pph(unittest.TestCase):
    def setUp(self):
        self.__input = os.path.join('data', 'snpeff', 'snpeff.vcf')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_snpeff2pph(self):
        sys.argv = ['', 'snpeff2bed', self.__input, self.__output]
        bioformats.cli.bioformats()
        sys.argv = ['', 'snpeff2bed', self.__input, self.__output,
                    '--bed3']
        bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
