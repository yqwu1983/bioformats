#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import bioformats.autosql
import os
import tempfile
import unittest
from bioformats.autosql import Table, TableEntry, Reader, Writer

try:
    import itertools.izip as zip
except ImportError:
    pass

path = os.path.dirname(__file__)
os.chdir(path)


class TestAutoSqlReader(unittest.TestCase):
    def setUp(self):
        self.__file = os.path.join(
            'data', 'autosql', 'correct.as'
        )

    def test_entries(self):
        """
        Check if the parser reads entries from an autoSql file.
        """
        with open(self.__file) as correct_file:
            reader = Reader(correct_file)
            for entry in reader.entries():
                self.assertIsInstance(entry, TableEntry)

    def test_get_table(self):
        """
        Check if the parser reads the whole table from an autoSql file.
        """
        with open(self.__file) as correct_file:
            reader = Reader(correct_file)
            self.assertIsInstance(reader.get_table(), Table)


class TestAutoSqlWriter(unittest.TestCase):
    def setUp(self):
        self.__file = os.path.join(
            'data', 'autosql', 'correct.as'
        )
        self.__output = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Test if an autoSql file is written correctly.
        """
        with open(self.__file) as input_file:
            autosql_table = Reader(input_file).get_table()

        with Writer(self.__output, autosql_table.name,
                    autosql_table.desc) as autosql_writer:
            for i in autosql_table.entries:
                autosql_writer.write(i)

        # compare the written file to the original one
        with open(self.__file) as original_file:
            with open(self.__output) as written_file:
                for line1, line2 in zip(original_file, written_file):
                    self.assertEqual(line1.rstrip(), line2.rstrip())

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)


class TestAutoSqlTypeRoutines(unittest.TestCase):
    def test_is_int(self):
        self.assertEqual(bioformats.autosql.is_int('5'), True)
        self.assertEqual(bioformats.autosql.is_int('5.4'), False)
        self.assertEqual(bioformats.autosql.is_int('5fa'), False)

    def test_is_float(self):
        self.assertEqual(bioformats.autosql.is_float('5'), True)
        self.assertEqual(bioformats.autosql.is_float('5.4'), True)
        self.assertEqual(bioformats.autosql.is_float('6.45sf'), False)

    def test_get_int_type(self):
        correct_types = [
            (255, 'ubyte'),
            (-128, 'byte'),
            (-129, 'short'),
            (256, 'ushort'),
            (65535, 'ushort'),
            (-32768, 'short'),
            (65537, 'uint'),
            (-32769, 'int'),
            (4294967295, 'uint'),
            (-2147483648, 'int')
        ]

        for x, y in correct_types:
            self.assertEqual(bioformats.autosql.get_int_type(x), y)

        self.assertIsNone(bioformats.autosql.get_int_type(pow(2, 32)))
        self.assertIsNone(bioformats.autosql.get_int_type(
            -pow(2, 31) - 1))

    def test_get_autosql_type(self):
        correct_types = [
            ('255', 'ubyte'),
            ('1.5', 'float'),
            ('a' * 255, 'string'),
            ('a' * 256, 'lstring')
        ]

        for x, y in correct_types:
            self.assertEqual(bioformats.autosql.get_autosql_type(x), y)

    def test_compare_autosql_types(self):
        triplets = [
            ('ubyte', 'byte', 'short'),
            ('uint', 'byte', 'int'),
            ('int', 'float', 'float'),
            ('string', 'byte', 'string'),
            ('lstring', 'string', 'lstring')
        ]

        for x in triplets:
            self.assertEqual(
                bioformats.autosql.compare_autosql_types(x[0], x[1]),
                x[2])
