#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import os
import logging
import tempfile
import unittest
from bioformats.gff3 import Record, Reader, Writer
from bioformats.exception import Gff3Error
try:
    import itertools.izip as zip
except ImportError:
    pass

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
        self.__incorrect_files = [
            os.path.join(self.__incorrect_file_dir, x) for x in
            glob.glob1(self.__incorrect_file_dir, '*.gff')]
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_records(self):
        """
        Check if the parse reads a file in the GFF3 format in the
        correct way.
        """
        # test against the correct input file
        with open(self.__correct_file) as gff_file:
            parser = Reader(gff_file)
            for record in parser.records():
                self.assertIsInstance(record, Record)
        # test against incorrect input files
        for i in self.__incorrect_files:
            with open(i) as gff_file:
                parser = Reader(gff_file)
                with self.assertRaises(Gff3Error):
                    for _ in parser.records():
                        pass


class TestGff3Writer(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join('data', 'gff3', 'correct.gff')
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Check if GFF3 lines are written properly to the output file.
        """
        with open(self.__input_file) as input_file:
            test_input = Reader(input_file)
            with Writer(self.__output_file) as test_output:
                for record in test_input.records():
                    test_output.write(record)

            # compare the test output file to the original one
            with open(self.__input_file) as test_input:
                with open(self.__output_file) as test_output:
                    for input_line, output_line in zip(test_input,
                                                       test_output):
                        self.assertEqual(input_line, output_line)
