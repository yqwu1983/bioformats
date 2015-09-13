#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
import bioformats.bed
from bioformats.repeatmasker import Reader, Record, whiten_color, \
    rmout2bed_record
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


class TestWhitenColor(unittest.TestCase):
    def test_whiten_color(self):
        self.assertEqual(whiten_color((0, 0, 0), 0), (0, 0, 0))
        self.assertEqual(whiten_color((0, 0, 0), 100), (255, 255, 255))
        self.assertEqual(whiten_color((50, 50, 50), 0), (50, 50, 50))
        self.assertEqual(whiten_color((255, 255, 255), 0), (255, 255,
                                                            255))
        self.assertEqual(whiten_color((255, 255, 255), 100), (255,
                                                              255, 255))
        with self.assertRaises(RepeatMaskerError):
            whiten_color((0, 0, 0), 255)


class TestRmOut2BedRecord(unittest.TestCase):
    def setUp(self):
        correct_file = os.path.join(
            'data', 'repeatmasker', 'repeatmasker.out'
        )
        with open(correct_file) as out_file:
            reader = Reader(out_file)
            self.__record = next(reader.repeats())

    def test_rmout2bed_record(self):
        self.assertIsInstance(rmout2bed_record(self.__record),
                              bioformats.bed.Record)
        # test the short option
        self.assertIsInstance(rmout2bed_record(self.__record,
                                               is_short=True),
                              bioformats.bed.Record)
        # test various values of the color argument
        for color in ('class', 'identity', 'class_identity'):
            self.assertIsInstance(
                rmout2bed_record(self.__record, color=color),
                bioformats.bed.Record
            )
        # test various values of the name argument
        for name in ('id', 'name', 'class', 'family', 'class_family'):
            self.assertIsInstance(
                rmout2bed_record(self.__record, name=name),
                bioformats.bed.Record
            )
        # test against incorrect values of the name and color arguments
        with self.assertRaises(RepeatMaskerError):
            rmout2bed_record(self.__record, name='error!')
        with self.assertRaises(RepeatMaskerError):
            rmout2bed_record(self.__record, color='error!')
        # test against an RNA-based repeat class
        rmout2bed_record(
            self.__record._replace(repeat_class_family='snRNA'),
            name='class'
        )
        # test against a missing repeat class for the option name=class
        rmout2bed_record(
            self.__record._replace(repeat_class_family='ERROR'),
            name='class'
        )
