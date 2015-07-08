#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import itertools
import logging
import os
import tempfile
import unittest
import bioformats.seqname
from bioformats.exception import IncorrectDictError
from bioformats.exception import MissingSeqNameError
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
        self.assertEqual(renamer.renaming_dict,
                         produced_renamer.renaming_dict)

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
                output_fasta.write(line)

        # perform the reverse renaming
        rev_renamer = bioformats.seqname.FastaSeqRenamer()
        rev_renamer.read_renaming_dict(self.__renaming_dict)
        with open(self.__rev_output, 'w') as rev_output_fasta:
            for line in renamer.renamed(self.__output, reverse=True):
                rev_output_fasta.write(line)

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


class TestNcbiFastaSeqRenamer(unittest.TestCase):
    def setUp(self):
        self.__output = tempfile.NamedTemporaryFile().name
        self.__test_dir = os.path.join(
            'data', 'seqrename',
        )
        self.__chr = os.path.join(
            self.__test_dir, 'chr_accessions_GRCh38.p2'
        )
        self.__unlocalized = os.path.join(
            self.__test_dir, 'unlocalized_accessions_GRCh38.p2'
        )
        self.__unplaced = os.path.join(
            self.__test_dir, 'unplaced_accessions_GRCh38.p2'
        )
        self.__acc_num_files = (self.__chr, self.__unlocalized,
                                self.__unplaced)

    def test_renamed(self):
        formats = ('refseq_full', 'genbank_full', 'refseq_gi',
                   'genbank_gi', 'refseq_acc_num', 'genbank_acc_num')
        for i, j in itertools.product(formats, repeat=2):
            renamer = bioformats.seqname.NcbiFastaSeqRenamer()
            for k in self.__acc_num_files:
                renamer.read_ncbi_acc_num(k, i, j)
            # convert sequence IDs
            input_file = os.path.join(self.__test_dir,
                                      'ncbi_' + i + '.fa')
            with open(self.__output, 'w') as output_file:
                for line in renamer.renamed(input_file):
                    output_file.write(line)

            example_file = os.path.join(self.__test_dir,
                                        'ncbi_' + j + '.fa')
            output_fasta = Fasta(self.__output)
            example_fasta = Fasta(example_file)
            # compare the obtained file to the example
            self.assertEqual(output_fasta.keys, example_fasta.keys)

        # remove temporary files and FASTA indices
        os.unlink(self.__output)
        os.unlink(self.__output + '.fai')
        for i in formats:
            os.unlink(os.path.join(self.__test_dir,
                                   'ncbi_' + i + '.fa.fai'))
