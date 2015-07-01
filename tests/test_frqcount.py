#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import logging
import os
import tempfile
import unittest
from bioformats.vcftools.frqcount import FrqCountRecord
from bioformats.vcftools.frqcount import Reader
from bioformats.vcftools.frqcount import SortedReader
from bioformats.vcftools.frqcount import Writer
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
        self.__compressed_correct_file = os.path.join(
            'data', 'frqcount', 'correct.gz'
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
        # test against the compressed correct input file
        parser = Reader(self.__compressed_correct_file, gzipped=True)
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
        for gzipped in (True, False):
            test_input = Reader(self.__input_file)
            with Writer(self.__output_file, gzipped) as test_output:
                for record in test_input.variants():
                    test_output.write(record)

            # compare the test output file to the original one
            test_output = Reader(self.__output_file, gzipped)
            for x, y in zip(
                    test_input.variants(), test_output.variants()):
                self.assertEqual(x, y)

            os.unlink(self.__output_file)


class TestSortedReader(unittest.TestCase):
    def setUp(self):
        self.__sorted = os.path.join(
            'data', 'frqcount', 'sorted.txt'
        )
        self.__unsorted = os.path.join(
            'data', 'frqcount', 'unsorted.txt'
        )
        self.__output = tempfile.NamedTemporaryFile().name

    def test_variants(self):
        """
        Check if SortedReader reads and sorts variants.
        """
        reader = SortedReader(self.__unsorted)
        with Writer(self.__output) as sorted_output:
            for chromosome in reader.variants:
                for variant in reader.variants[chromosome]:
                    sorted_output.write(variant)

        # compare the original sorted file and the produced one
        test_original = Reader(self.__sorted)
        test_output = Reader(self.__output)
        for x, y in zip(
                test_original.variants(), test_output.variants()):
            self.assertEqual(x, y)

        os.unlink(self.__output)
