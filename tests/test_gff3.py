#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import os
import logging
import tempfile
import unittest
from bioformats.gff3 import Record, Reader, Writer, analyze_tags
from bioformats.gff3 import gff2to3
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

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)


class TestAnalyzeTags(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join('data', 'gff3', 'correct.gff')

    def test_analyze_tags(self):
        """
        Check if tags are analyzed in a proper way.
        """
        with open(self.__input_file) as gff_file:
            result = analyze_tags(gff_file)
            self.assertIsInstance(result, dict)

        with open(self.__input_file) as gff_file:
            result = analyze_tags(gff_file, feature_source='example')
            self.assertIsInstance(result, dict)

        with open(self.__input_file) as gff_file:
            result = analyze_tags(gff_file, feature_type='CDS')
            self.assertIsInstance(result, dict)


class TestGff2to3(unittest.TestCase):
    def setUp(self):
        self.__gff2 = os.path.join('data', 'gff3', 'example.gff2')
        self.__incorrect_gff2 = os.path.join('data', 'gff3',
                                             'incorrect.gff2')
        self.__gff3 = os.path.join('data', 'gff3', 'example.gff3')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_gff2to3(self):
        with open(self.__gff2) as input_file:
            with open(self.__output, 'w') as output_file:
                gff2to3(input_file, output_file)

        # compare the produced GFF3 file to the correct one
        with open(self.__output) as produced_gff:
            with open(self.__gff3) as correct_gff:
                for x, y in zip(produced_gff, correct_gff):
                    self.assertEqual(x.rstrip(), y.rstrip())

        # try the incorrect GFF2 file in the strict mode
        with open(self.__incorrect_gff2) as input_file:
            with open(self.__output, 'w') as output_file:
                with self.assertRaises(Gff3Error):
                    gff2to3(input_file, output_file, strict=True)

        # try the incorrect GFF2 file ignoring incorrect input lines
        with open(self.__incorrect_gff2) as input_file:
            with open(self.__output, 'w') as output_file:
                gff2to3(input_file, output_file, strict=False)

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
