#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli
from bioformats.seqname import BaseSeqRenamer
from pyfaidx import Fasta
try:
    import itertools.izip as zip
except ImportError:
    pass


class TestNcbiRenameSeq(unittest.TestCase):
    def setUp(self):
        self.__table = os.path.join(
            'data', 'seqrename', 'ncbi_table.txt'
        )
        self.__fasta = os.path.join(
            'data', 'seqrename', 'ncbi_genbank.fa'
        )
        self.__ucsc_fasta = os.path.join(
            'data', 'seqrename', 'ncbi_ucsc.fa'
        )
        self.__chr = os.path.join(
            'data', 'seqrename', 'chr_accessions_GRCh38.p2'
        )
        self.__unloc = os.path.join(
            'data', 'seqrename', 'unlocalized_accessions_GRCh38.p2'
        )
        self.__unpl = os.path.join(
            'data', 'seqrename', 'unplaced_accessions_GRCh38.p2'
        )
        self.__output_table = tempfile.NamedTemporaryFile().name
        self.__output = tempfile.NamedTemporaryFile().name
        self.__rev_output = tempfile.NamedTemporaryFile().name

    def test_ncbiseqrename_fasta(self):
        """
        Check if NCBI sequence names in a FASTA file are properly
        changed.
        """
        sys.argv = ['', self.__fasta, 'genbank', self.__output, 'ucsc',
                    '--chr', self.__chr, '--unloc', self.__unloc,
                    '--unpl', self.__unpl, '--fasta']

        bioformats.cli.ncbirenameseq()

        # check if the obtained and original files are the same
        original_fasta = Fasta(self.__ucsc_fasta)
        renamed_fasta = Fasta(self.__output)
        for x, y in zip(original_fasta.keys(), renamed_fasta.keys()):
            self.assertEqual(x, y)

        os.unlink(self.__ucsc_fasta + '.fai')
        os.unlink(self.__output)
        os.unlink(self.__output + '.fai')

    def test_ncbiseqrename_table(self):
        """
        Check if NCBI sequence names in a table are properly changed.
        """
        sys.argv = ['', self.__table, 'genbank', self.__output,
                    'refseq',
                    '--chr', self.__chr, '--unloc', self.__unloc,
                    '--unpl', self.__unpl,
                    '--output_table', self.__output_table]

        bioformats.cli.ncbirenameseq()

        # check if the produced renaming table is correct
        test_renamer = BaseSeqRenamer()
        test_renamer.read_renaming_dict(self.__output_table)

        sys.argv = ['', self.__output, 'refseq', self.__rev_output,
                    'genbank',
                    '--chr', self.__chr, '--unloc', self.__unloc,
                    '--unpl', self.__unpl]

        bioformats.cli.ncbirenameseq()

        # check if the obtained and original files are the same
        # compare the original and reverse-renamed tables
        with open(self.__table) as original_table:
            with open(self.__rev_output) as rev_renamed_table:
                for x, y in zip(original_table, rev_renamed_table):
                    self.assertEqual(x, y)

        os.unlink(self.__output_table)
        os.unlink(self.__output)
        os.unlink(self.__rev_output)
