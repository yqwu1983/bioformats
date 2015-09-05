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
        correct_bed = Reader(self.__correct_bed)
        output_bed = Reader(self.__output_file)
        for x, y in zip(correct_bed.records(), output_bed.records()):
            self.assertEqual(x, y)

        os.unlink(self.__output_file)
        os.unlink(self.__fasta + '.fai')
