#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import os
import tempfile
import unittest
import bioformats.seqname
from bioformats.exception import IncorrectDictError
from bioformats.exception import MissingSeqNameError
from future.utils import itervalues
from pyfaidx import Fasta
try:
    import itertools.izip as zip
except ImportError:
    pass


class TestBaseSeqRenamer(unittest.TestCase):
    def setUp(self):
        self.__correct = os.path.join(
            'data', 'seqrename', 'correct_dict.txt'
        )
        self.__incorrect = os.path.join(
            'data', 'seqrename', 'incorrect_dict.txt'
        )
        self.__output = tempfile.NamedTemporaryFile().name

        # silence logging messages
        logging.disable(logging.ERROR)

    def test_read(self):
        """
        Check if a renaming dictionary is correctly read.
        """
        # check against a correct input file
        renamer = bioformats.seqname.BaseSeqRenamer()
        renamer.read_renaming_dict(self.__correct)
        # check against an incorrect input file
        with self.assertRaises(IncorrectDictError):
            renamer.read_renaming_dict(self.__incorrect)

    def test_write(self):
        """
        Check if a renaming dictionary is correctly written.
        """
        renamer = bioformats.seqname.BaseSeqRenamer()
        renamer.read_renaming_dict(self.__correct)
        renamer.write_renaming_dict(self.__output)
        # compare the original and written dictionaries
        produced_renamer = bioformats.seqname.BaseSeqRenamer()
        produced_renamer.read_renaming_dict(self.__output)
        for x, y in zip(
                itervalues(renamer.renaming_dict),
                itervalues(produced_renamer.renaming_dict)):
            self.assertEqual(x, y)

        os.unlink(self.__output)


class TestFastaSeqRenamer(unittest.TestCase):
    def setUp(self):
        self.__renaming_dict = os.path.join(
            'data', 'seqrename', 'correct_dict.txt'
        )
        self.__fasta = os.path.join(
            'data', 'seqrename', 'fasta.fa'
        )
        self.__output = tempfile.NamedTemporaryFile().name
        self.__rev_output = tempfile.NamedTemporaryFile().name

    def tearDown(self):
        # remove the FASTA index file created while we iterated
        # through its sequences
        os.unlink(self.__fasta + '.fai')

    def test_renamed(self):
        """
        Check if sequences in a FASTA file are properly renamed.
        """
        renamer = bioformats.seqname.FastaSeqRenamer()
        renamer.read_renaming_dict(self.__renaming_dict)
        with open(self.__output, 'w') as output_fasta:
            for line in renamer.renamed(self.__fasta):
                output_fasta.write(line + '\n')

        # perform the reverse renaming
        rev_renamer = bioformats.seqname.FastaSeqRenamer()
        rev_renamer.read_renaming_dict(self.__renaming_dict)
        with open(self.__rev_output, 'w') as rev_output_fasta:
            for line in renamer.renamed(self.__output, reverse=True):
                rev_output_fasta.write(line + '\n')

        # compare the original and reverse-renamed FASTA files
        original_fasta = Fasta(self.__fasta)
        rev_renamed_fasta = Fasta(self.__rev_output)
        for x, y in zip(
                original_fasta.keys(),
                rev_renamed_fasta.keys()):
            self.assertEqual(x, y)

        # check if the missing sequence exception is raised
        del renamer.renaming_dict['seq2']
        with self.assertRaises(MissingSeqNameError):
            for _ in renamer.renamed(self.__fasta):
                pass

        os.unlink(self.__output)
        os.unlink(self.__rev_output)


class TestTableSeqRenamer(unittest.TestCase):
    def setUp(self):
        self.__renaming_dict = os.path.join(
            'data', 'seqrename', 'correct_dict.txt'
        )
        self.__table = os.path.join(
            'data', 'seqrename', 'table.txt'
        )
        self.__output = tempfile.NamedTemporaryFile().name
        self.__rev_output = tempfile.NamedTemporaryFile().name

    def test_renamed(self):
        renamer = bioformats.seqname.TableSeqRenamer()
        renamer.read_renaming_dict(self.__renaming_dict)
        with open(self.__output, 'w') as output_table:
            for line in renamer.renamed(self.__table, 0):
                output_table.write(line + '\n')

        # perform the reverse renaming
        rev_renamer = bioformats.seqname.TableSeqRenamer()
        rev_renamer.read_renaming_dict(self.__renaming_dict)
        with open(self.__rev_output, 'w') as rev_output_table:
            for line in renamer.renamed(self.__output, 0, reverse=True):
                rev_output_table.write(line + '\n')

        # compare the original and reverse-renamed tables
        with open(self.__table) as original_table:
            with open(self.__rev_output) as rev_renamed_table:
                for x, y in zip(original_table, rev_renamed_table):
                    self.assertEqual(x, y)

        # check if the missing sequence exception is raised
        del renamer.renaming_dict['seq2']
        with self.assertRaises(MissingSeqNameError):
            for _ in renamer.renamed(self.__table, 0):
                pass

        os.unlink(self.__output)
        os.unlink(self.__rev_output)
