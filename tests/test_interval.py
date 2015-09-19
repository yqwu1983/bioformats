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
        parser = Reader(self.__correct_file)
        for record in parser.intervals():
            self.assertIsInstance(record, Record)
        # test against the incorrect input file
        parser = Reader(self.__incorrect_file)
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
        test_input = Reader(self.__input_file)
        with Writer(self.__output_file) as test_output:
            for record in test_input.intervals():
                test_output.write(record)

        # compare the test output file to the original one
        test_output = Reader(self.__output_file)
        for x, y in zip(
                test_input.intervals(), test_output.intervals()):
            self.assertEqual(x, y)

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)
