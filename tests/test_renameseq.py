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

path = os.path.dirname(__file__)
os.chdir(path)


class TestRenameSeq(unittest.TestCase):
    def setUp(self):
        self.__renaming_dict = os.path.join(
            'data', 'seqrename', 'correct_dict.txt'
        )
        self.__fasta = os.path.join(
            'data', 'seqrename', 'fasta.fa'
        )
        self.__fasta_desc = os.path.join(
            'data', 'seqrename', 'fasta_desc.fa'
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
        sys.argv = ['', 'renameseq', self.__renaming_dict,
                    self.__fasta, self.__output, '-f']

        bioformats.cli.bioformats()

        sys.argv = ['', 'renameseq', self.__renaming_dict,
                    self.__output, self.__rev_output, '-f', '-r']

        bioformats.cli.bioformats()

        # check if the obtained and original files are the same
        original_fasta = Fasta(self.__fasta)
        rev_renamed_fasta = Fasta(self.__rev_output)
        for x, y in zip(
                original_fasta.keys(),
                rev_renamed_fasta.keys()):
            self.assertEqual(x, y)

        # check if sequence descriptions are removed
        sys.argv = ['', 'renameseq', self.__renaming_dict,
                    self.__fasta, self.__rev_output, '-f',
                    '--no_description']

        bioformats.cli.bioformats()

        with open(self.__output) as renamed_fasta:
            with open(self.__rev_output) as nodesc_renamed_fasta:
                for x, y in zip(renamed_fasta, nodesc_renamed_fasta):
                    self.assertEqual(x, y)

    def test_renameseq_table(self):
        """
        Check if sequence names in a tabular file are properly changed.
        """
        sys.argv = ['', 'renameseq', self.__renaming_dict,
                    self.__table, self.__output]

        bioformats.cli.bioformats()

        sys.argv = ['', 'renameseq', self.__renaming_dict,
                    self.__output, self.__rev_output, '-r']

        bioformats.cli.bioformats()

        # check if the obtained and original files are the same
        # compare the original and reverse-renamed tables
        with open(self.__table) as original_table:
            with open(self.__rev_output) as rev_renamed_table:
                for x, y in zip(original_table, rev_renamed_table):
                    self.assertEqual(x, y)

    def tearDown(self):
        for i in (self.__output, self.__rev_output,
                  self.__output + '.fai', self.__rev_output + '.fai'):
            if os.path.isfile(i):
                os.unlink(i)
