#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli
from pyfaidx import Fasta
try:
    import itertools.izip as zip
except ImportError:
    pass


class TestRenameSeq(unittest.TestCase):
    def setUp(self):
        self.__renaming_dict = os.path.join(
            'data', 'seqrename', 'correct_dict.txt'
        )
        self.__fasta = os.path.join(
            'data', 'seqrename', 'fasta.fa'
        )
        self.__table = os.path.join(
            'data', 'seqrename', 'table.txt'
        )
        self.__output = tempfile.NamedTemporaryFile().name
        self.__rev_output = tempfile.NamedTemporaryFile().name

    def test_renameseq_fasta(self):
        """
        Check if sequence names in a FASTA file are properly changed.
        """
        sys.argv = ['', self.__renaming_dict, self.__fasta,
                    self.__output, '-f']

        bioformats.cli.renameseq()

        sys.argv = ['', self.__renaming_dict, self.__output,
                    self.__rev_output, '-f', '-r']

        bioformats.cli.renameseq()

        # check if the obtained and original files are the same
        original_fasta = Fasta(self.__fasta)
        rev_renamed_fasta = Fasta(self.__rev_output)
        for x, y in zip(
                original_fasta.keys(),
                rev_renamed_fasta.keys()):
            self.assertEqual(x, y)

        os.unlink(self.__output)
        os.unlink(self.__rev_output)

    def test_renameseq_table(self):
        """
        Check if sequence names in a tabular file are properly changed.
        """
        sys.argv = ['', self.__renaming_dict, self.__table,
                    self.__output]

        bioformats.cli.renameseq()

        sys.argv = ['', self.__renaming_dict, self.__output,
                    self.__rev_output, '-r']

        bioformats.cli.renameseq()

        # check if the obtained and original files are the same
        # compare the original and reverse-renamed tables
        with open(self.__table) as original_table:
            with open(self.__rev_output) as rev_renamed_table:
                for x, y in zip(original_table, rev_renamed_table):
                    self.assertEqual(x, y)

        os.unlink(self.__output)
        os.unlink(self.__rev_output)
