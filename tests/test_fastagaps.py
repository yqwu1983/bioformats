#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli
from bioformats.bed import Reader
try:
    import itertools.izip as zip
except ImportError:
    pass


class TestFastaGaps(unittest.TestCase):
    def setUp(self):
        self.__fasta = os.path.join(
            'data', 'fastagaps', 'gaps.fa'
        )
        self.__correct_bed = os.path.join(
            'data', 'fastagaps', 'gaps.bed'
        )
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_fastagaps(self):
        """
        Test if gaps are correctly identified in a FASTA file.
        """
        sys.argv = ['', 'fastagaps', self.__fasta, self.__output_file]
        bioformats.cli.bioformats()

        # compare the obtained BED file with the correct one
        with open(self.__correct_bed) as correct_bed:
            with open(self.__output_file) as output_bed:
                correct_reader = Reader(correct_bed)
                output_reader = Reader(output_bed)
                for x, y in zip(
                        correct_reader.records(),
                        output_reader.records()):
                    self.assertEqual(x, y)

    def tearDown(self):
        for i in (self.__output_file, self.__fasta + '.fai'):
            if os.path.isfile(i):
                os.unlink(i)
