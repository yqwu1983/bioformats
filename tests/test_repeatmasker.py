#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from bioformats.repeatmasker import Reader, Record
from bioformats.exception import RepeatMaskerError


class TestRepeatMaskerReader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'repeatmasker', 'repeatmasker.out'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'repeatmasker', 'incorrect'
        )
        self.__incorrect_files = [os.path.join(
            self.__incorrect_file_dir, x) for x in os.listdir(
            self.__incorrect_file_dir)]
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_records(self):
        # test againt the correct input file
        with open(self.__correct_file) as correct_file:
            reader = Reader(correct_file)
            for repeat in reader.repeats():
                self.assertIsInstance(repeat, Record)
        # test against incorrect input files
        for i in self.__incorrect_files:
            with open(i) as incorrect_file:
                    reader = Reader(incorrect_file)
                    with self.assertRaises(RepeatMaskerError):
                        for _ in reader.repeats():
                            pass
