#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import tempfile
import unittest
import bioformats.autosql
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
        self.__correct_file_names = (
            'correct.bed12',
            'correct_aux.bed6'
        )
        self.__column_test_file = os.path.join(
            'data', 'bed', 'column_test.bed'
        )
        self.__correct_columns = (
            (12, 1),
            (9, 4),
            (9, 4),
            (9, 4),
            (8, 5),
            (6, 7),
            (6, 7),
            (5, 8),
            (4, 9),
        )
        self.__correct_files = [os.path.join('data', 'bed', x)
                                for x in self.__correct_file_names]
        self.__incorrect_file = os.path.join(
            'data', 'bed', 'incorrect_aux.bed6'
        )
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_records(self):
        """
        Check if the parser reads a file in the BED format in the
        correct way.
        """
        # test against the correct input files
        for i in self.__correct_files:
            with open(i) as bed_file:
                parser = Reader(bed_file)
                for record in parser.records():
                    self.assertIsInstance(record, Record)
        # test against an unsorted input file
        for i in self.__correct_files:
            with open(i) as unsorted_bed_file:
                parser = Reader(unsorted_bed_file)
                with self.assertRaises(BedError):
                    for record in parser.records(check_order=True):
                        self.assertIsInstance(record, Record)
        # test against the incorrect input file
        with open(self.__incorrect_file) as bed_file:
            parser = Reader(bed_file)
            with self.assertRaises(BedError):
                for record in parser.records():
                    self.assertIsInstance(record, Record)

    def test_columns(self):
        """
        Check if BED and auxiliary columns are correctly counted.
        """
        with open(self.__column_test_file) as bed_file:
            reader = Reader(bed_file)
            for i, record in enumerate(reader.records()):
                self.assertEqual(reader.bed_columns,
                                 self.__correct_columns[i][0])
                self.assertEqual(reader.aux_columns,
                                 self.__correct_columns[i][1])


class TestBedWriter(unittest.TestCase):
    def setUp(self):
        self.__input_file_names = (
            'correct.bed12',
            'correct_aux.bed6'
        )
        self.__input_file = [os.path.join('data', 'bed', x)
                             for x in self.__input_file_names]
        self.__output_file = tempfile.NamedTemporaryFile().name
        
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_write(self):
        """
        Check if BED records are written in the correct way.
        """
        for i in self.__input_file:
            with open(i) as bed_file:
                bed_input = Reader(bed_file)
                with Writer(self.__output_file) as bed_output:
                    for record in bed_input.records():
                        bed_output.write(record)

            # check if the lines are identical
            with open(i) as original_file, \
                    open(self.__output_file) as written_file:
                for x, y in zip(original_file, written_file):
                    self.assertEqual(x.rstrip(), y.rstrip())

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


class TestGetAutoSqlTable(unittest.TestCase):
    def setUp(self):
        self.__input_file_names = (
            'correct.bed12',
            'correct_aux.bed6'
        )
        self.__input_file = [os.path.join('data', 'bed', x)
                             for x in self.__input_file_names]

    def test_get_autosql_table(self):
        for i in self.__input_file:
            with open(i) as bed_file:
                reader = Reader(bed_file)
                autosql_table = bioformats.bed.get_autosql_table(
                    reader
                )
                self.assertIsInstance(autosql_table,
                                      bioformats.autosql.Table)
                self.assertEqual(len(autosql_table.entries),
                                 reader.bed_columns +
                                 reader.aux_columns)


class TestGetBlocks(unittest.TestCase):
    def test_get_blocks(self):
        with self.assertRaises(BedError):
            bioformats.bed.get_blocks([1, 5], [3])
        start = [1, 5, 18]
        end = [3, 10, 29]
        result = bioformats.bed.get_blocks(start, end)
        self.assertEqual(result[0], [0, 4, 17])
        self.assertEqual(result[1], [3, 6, 12])


class TestConvertGff2BedGene(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join(
            'data', 'gff3', 'exons.gff'
        )
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_convert_gff2bed_gene(self):
        bioformats.bed.convert_gff2bed_gene(self.__input_file,
                                            self.__output_file)
        # try to read the obtained BED file
        with open(self.__output_file) as bed_file:
            reader = Reader(bed_file)
            for _ in reader.records():
                pass

    def test_convert_gff2bed(self):
        bioformats.bed.convert_gff2bed(self.__input_file,
                                       self.__output_file,
                                       'exon')
        # try to read the obtained BED file
        with open(self.__output_file) as bed_file:
            reader = Reader(bed_file)
            for _ in reader.records():
                pass
        bioformats.bed.convert_gff2bed(self.__input_file,
                                       self.__output_file,
                                       'exon',
                                       attributes=['Parent'],
                                       name_tag='ID')
        # try to read the obtained BED file
        with open(self.__output_file) as bed_file:
            reader = Reader(bed_file)
            for _ in reader.records():
                pass
        bioformats.bed.convert_gff2bed(self.__input_file,
                                       self.__output_file,
                                       'exon',
                                       attributes=['Parent'],
                                       name_tag='NoTag')
        # try to read the obtained BED file
        with open(self.__output_file) as bed_file:
            reader = Reader(bed_file)
            for _ in reader.records():
                pass

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)
