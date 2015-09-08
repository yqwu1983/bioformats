#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli


class TestBedAutoSql(unittest.TestCase):
    def setUp(self):
        self.__bed_file_names = (
            'correct_aux.bed6',
            'correct.bed12',
            'incorrect_aux.bed6'
        )
        self.__bed_files = [os.path.join('data', 'bed', x)
                            for x in self.__bed_file_names]
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_bedcolumns(self):
        """
        Test if the bedcolumns tool correctly processes BED files.
        """
        for i in self.__bed_files:
            sys.argv = ['', 'bedautosql', i, self.__output_file]
            bioformats.cli.bioformats()
