#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

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