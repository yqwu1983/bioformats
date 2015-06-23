#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from bioformats.bed import BedRecord, Reader, Writer
from bioformats.exception import BedError

try:
    import itertools.izip as zip
except ImportError:
    pass

path = os.path.dirname(__file__)
os.chdir(path)


class TestBedReader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'bed', 'correct.bed'
        )
        self.__incorrect_file = os.path.join(
            'data', 'bed', 'incorrect.bed'
        )
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_records(self):
        """
        Check if the parser reads a file in the BED format in the
        correct way.
        """
        # test against the correct input file
        parser = Reader(self.__correct_file)
        for record in parser.records():
            self.assertIsInstance(record, BedRecord)
        # test against the incorrect input file
        parser = Reader(self.__incorrect_file)
        with self.assertRaises(BedError):
            for record in parser.records():
                self.assertIsInstance(record, BedRecord)


class TestBedWriter(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join(
            'data', 'bed', 'correct.bed'
        )
        self.__output_file = os.path.join(
            'data', 'bed', 'test.bed'
        )
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_write(self):
        """
        Check if BED records are written in the correct way.
        """
        bed_input = Reader(self.__input_file)
        with Writer(self.__output_file) as bed_output:
            for record in bed_input.records():
                bed_output.write(record)

        # check if the lines are identical
        with open(self.__input_file) as original_file, \
                open(self.__output_file) as written_file:
            for x, y in zip(original_file, written_file):
                self.assertEqual(x, y)

        os.unlink(self.__output_file)
