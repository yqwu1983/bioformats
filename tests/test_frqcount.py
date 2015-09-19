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
        with open(self.__correct_file) as input_file:
            for record in Reader(input_file).variants():
                self.assertIsInstance(record, FrqCountRecord)
        # test against incorrect input files
        for frqcount_file in self.__incorrect_files:
            with open(os.path.join(self.__incorrect_file_dir,
                                   frqcount_file)) as input_file:
                with self.assertRaises(FrqCountReaderError):
                    for _ in Reader(input_file).variants():
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
        with open(self.__input_file) as input_file, \
                open(self.__output_file, 'w') as output_file:
            writer = Writer(output_file)
            for record in Reader(input_file).variants():
                writer.write(record)

        # compare the test output file to the original one
        with open(self.__input_file) as input_file, \
                open(self.__output_file) as output_file:
            for x, y in zip(
                    Reader(input_file).variants(),
                    Reader(output_file).variants()):
                self.assertEqual(x, y)

    def tearDown(self):
        if os.path.isfile(self.__output_file):
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
        with open(self.__unsorted) as unsorted_input, \
                open(self.__output, 'w') as sorted_output:
            reader = SortedReader(unsorted_input)
            writer = Writer(sorted_output)
            for chromosome in reader.variants:
                for variant in reader.variants[chromosome]:
                    writer.write(variant)

        # compare the original sorted file and the produced one
        with open(self.__sorted) as original_sorted, \
                open(self.__output) as produced_sorted:
            for x, y in zip(
                    Reader(original_sorted).variants(),
                    Reader(produced_sorted).variants()):
                self.assertEqual(x, y)

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
