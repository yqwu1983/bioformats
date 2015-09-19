#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import os
import logging
import unittest
from bioformats.lav import Lav
from bioformats.exception import LavError


class TestLavAlignment(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'lav', 'lav_alignments.txt'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'lav', 'incorrect_input'
        )
        self.__incorrect_files = glob.glob1(self.__incorrect_file_dir,
                                            '*.txt')
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_alignments(self):
        """
        Check if the parser reads a file in the LAV format in the
        correct way.
        """
        # test against the correct input file
        parser = Lav(self.__correct_file)
        for alignment in parser.alignments():
            self.assertEqual(len(alignment), 7)
        for alignment in parser.alignments(gapped=False):
            self.assertEqual(len(alignment), 8)
        # test againts incorrect input files
        for lav_file in self.__incorrect_files:
            parser = Lav(os.path.join(self.__incorrect_file_dir,
                                      lav_file))
            with self.assertRaises(LavError):
                for alignment in parser.alignments():
                    self.assertIsInstance(alignment,
                                          Lav.GapFreeAlignment)
