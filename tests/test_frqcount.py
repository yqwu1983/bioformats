#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import logging
import os
import tempfile
import unittest
from bioformats.vcftools.frqcount import Reader, Writer, FrqCountRecord
from bioformats.exception import FrqCountReaderError
try:
    import itertools.izip as zip
except ImportError:
    pass

path = os.path.dirname(__file__)
os.chdir(path)


class TestReader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'frqcount', 'correct.txt'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'frqcount', 'incorrect_input'
        )
        self.__incorrect_files = glob.glob1(
            self.__incorrect_file_dir, '*.txt')

        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_variants(self):
        """
        Check if the parser reads a file in the VCFtools frequency
        counts format in the correct way.
        """
        # test against the correct input file
        parser = Reader(self.__correct_file)
        for record in parser.variants():
            self.assertIsInstance(record, FrqCountRecord)
        # test against incorrect input files
        for frqcount_file in self.__incorrect_files:
            parser = Reader(os.path.join(self.__incorrect_file_dir,
                                         frqcount_file))
            with self.assertRaises(FrqCountReaderError):
                for _ in parser.variants():
                    pass


class TestWriter(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join('data', 'frqcount',
                                         'correct.txt')
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Check if VCFtools allele frequencies are correctly written to
        the output file.
        """
        test_input = Reader(self.__input_file)
        with Writer(self.__output_file) as test_output:
            for record in test_input.variants():
                test_output.write(record)

        # compare the test output file to the original one
        with open(self.__input_file) as test_input:
            with open(self.__output_file) as test_output:
                for input_line, output_line in zip(test_input,
                                                   test_output):
                    self.assertEqual(input_line, output_line)
