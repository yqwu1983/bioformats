#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import os
import tempfile
import unittest
from bioformats.interval import Record
from bioformats.interval import Reader
from bioformats.interval import Writer
from bioformats.exception import IntervalError
try:
    import itertools.izip as zip
except ImportError:
    pass

path = os.path.dirname(__file__)
os.chdir(path)


class TestReader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'interval', 'correct.txt'
        )
        self.__incorrect_file = os.path.join(
            'data', 'interval', 'incorrect.txt'
        )

        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_intervals(self):
        """
        Check if the parser reads a file in the interval format in
        the correct way.
        """
        # test against the correct input file
        with open(self.__correct_file) as input_file:
            parser = Reader(input_file)
            for record in parser.intervals():
                self.assertIsInstance(record, Record)
        # test against the incorrect input file
        with open(self.__incorrect_file) as input_file:
            parser = Reader(input_file)
            with self.assertRaises(IntervalError):
                for _ in parser.intervals():
                    pass


class TestWriter(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join('data', 'interval',
                                         'correct.txt')
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Check if intervals are correctly written to the output file.
        """
        with open(self.__input_file) as test_input:
            input_reader = Reader(test_input)
            with Writer(self.__output_file) as test_output:
                for record in input_reader.intervals():
                    test_output.write(record)

            # compare the test output file to the original one
        with open(self.__input_file) as test_input:
            with open(self.__output_file) as produced_output:
                input_reader = Reader(test_input)
                output_reader = Reader(produced_output)
                for x, y in zip(
                        input_reader.intervals(),
                        output_reader.intervals()):
                    self.assertEqual(x, y)

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)
