#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import tempfile
import unittest
import bioformats.bed
from bioformats.bed import Record, Reader, Writer
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
        with open(self.__correct_file) as bed_file:
            parser = Reader(bed_file)
            for record in parser.records():
                self.assertIsInstance(record, Record)
        # test against the incorrect input file
        with open(self.__incorrect_file) as bed_file:
            parser = Reader(bed_file)
            with self.assertRaises(BedError):
                for record in parser.records():
                    self.assertIsInstance(record, Record)


class TestBedWriter(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join(
            'data', 'bed', 'correct.bed'
        )
        self.__output_file = tempfile.NamedTemporaryFile().name
        
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_write(self):
        """
        Check if BED records are written in the correct way.
        """
        with open(self.__input_file) as bed_file:
            bed_input = Reader(bed_file)
            with Writer(self.__output_file) as bed_output:
                for record in bed_input.records():
                    bed_output.write(record)

        # check if the lines are identical
        with open(self.__input_file) as original_file, \
                open(self.__output_file) as written_file:
            for x, y in zip(original_file, written_file):
                self.assertEqual(x, y)

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)


class TestBedContentsCheckRoutines(unittest.TestCase):
    def test_is_score(self):
        test_values = [
            ('-1', False),
            ('0', True),
            ('1000', True),
            ('1001', False)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_score(x), y)

    def test_is_strand(self):
        test_values = [
            ('+', True),
            ('-', True),
            ('.', False)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_strand(x), y)

    def test_is_coord(self):
        test_values = [
            ('AB', False),
            ('-1', False),
            ('0', True),
            ('10', True)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_coord(x), y)

    def test_itemrgb(self):
        test_values = [
            ('60,35', False),
            ('255,0,B', False),
            ('255,-1,0', False),
            ('256,0,0', False),
            ('255,0,255', True)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_itemrgb(x), y)

    def test_is_block_count(self):
        test_values = [
            ('F', False),
            ('0', False),
            ('1', True)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_block_count(x), y)

    def test_is_block_sizes(self):
        test_values = [
            ('ABF', False),
            ('0,FG', False),
            ('1,0,2', False),
            ('1,5,2', True)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_block_sizes(x), y)

    def test_is_block_starts(self):
        test_values = [
            ('GHT', False),
            ('0', True),
            ('0,F', False),
            ('1,2', False),
            ('0,5,2', False),
            ('0,2', True)
        ]
        for x, y in test_values:
            self.assertEqual(bioformats.bed.is_block_starts(x), y)