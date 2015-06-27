#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import logging
import os
import unittest
from bioformats.vcftools.frqcount import Reader
from bioformats.exception import FrqCountReaderError

path = os.path.dirname(__file__)
os.chdir(path)


class TestReader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'frqcount', 'correct.txt'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'frqcount', 'incorrect_input'
        )
        self.__incorrect_files = glob.glob1(
            self.__incorrect_file_dir, '*.txt')

        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_variants(self):
        """
        Check if the parser reads a file in the VCFtools frequency
        counts format in the correct way.
        """
        # test against the correct input file
        parser = Reader(self.__correct_file)
        for record in parser.variants():
            self.assertIsInstance(record, Reader.Record)
        # test against incorrect input files
        for frqcount_file in self.__incorrect_files:
            parser = Reader(os.path.join(self.__incorrect_file_dir,
                                         frqcount_file))
            with self.assertRaises(FrqCountReaderError):
                for _ in parser.variants():
                    pass
