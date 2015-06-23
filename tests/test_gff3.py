#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import os
import logging
import tempfile
import unittest
from bioformats.gff3 import Gff3Record, Reader, Writer
from bioformats.exception import Gff3Error

path = os.path.dirname(__file__)
os.chdir(path)


class TestGff3Reader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'gff3', 'correct.gff'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'gff3', 'incorrect_input'
        )
        self.__incorrect_files = glob.glob1(
            self.__incorrect_file_dir, '*.gff')
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_records(self):
        """
        Check if the parse reads a file in the GFF3 format in the
        correct way.
        """
        # test against the correct input file
        parser = Reader(self.__correct_file)
        for record in parser.records():
            self.assertIsInstance(record, Gff3Record)
        # test against incorrect input files
        for gff_file in self.__incorrect_files:
            parser = Reader(os.path.join(self.__incorrect_file_dir,
                                         gff_file))
            with self.assertRaises(Gff3Error):
                for _ in parser.records():
                    pass


class TestGff3Writer(unittest.TestCase):
    def setUp(self):
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Check if GFF3 lines are written properly to the output file.
        """
        record1 = Gff3Record('ctg123', '.', 'exon', 1300, 1500,
                             100.0, '+', '.', dict(ID='exon00001'))
        record2 = Gff3Record('ctg123', '.', 'exon', 7000, 9000,
                             34.5, '+', '.', None)
        with Writer(self.__output_file) as output_gff:
            output_gff.write(record1)
            output_gff.write(record2)

        # read records from the file and compare it to the original
        # records
        reader = Reader(self.__output_file)
        written_records = [i for i in reader.records()]

        self.assertEqual(record1, written_records[0])
        self.assertEqual(record2, written_records[1])
