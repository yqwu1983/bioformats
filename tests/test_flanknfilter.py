#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestFlankNFilter(unittest.TestCase):
    def setUp(self):
        self.test_fa = os.path.join(
            'data', 'fasta', 'flanknfilter_test.fa'
        )
        self.test_bed = os.path.join(
            'data', 'fasta', 'flanknfilter_test.bed'
        )
        self.test_vcf = os.path.join(
            'data', 'fasta', 'flanknfilter_test.vcf'
        )
        self.output = tempfile.NamedTemporaryFile().name

    def test_flanknfilter(self):
        sys.argv = ['', 'flanknfilter', '-t', 'bed', self.test_bed,
                    self.test_fa, self.output]
        bioformats.cli.bioformats()
        sys.argv = ['', 'flanknfilter', '-t', 'vcf', self.test_vcf,
                    self.test_fa, self.output]
        bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.isfile(self.output):
            os.unlink(self.output)
