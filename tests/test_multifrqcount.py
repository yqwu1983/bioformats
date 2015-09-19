#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import os
import unittest
from bioformats.vcftools.multifrqcount import Reader
from bioformats.vcftools.multifrqcount import Record


class TestReader(unittest.TestCase):
    def setUp(self):
        self.__correct = os.path.join(
            'data', 'multifrqcount', 'correct.txt'
        )

        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_alleles(self):
        """
        Check if allele counts are read correctly.
        """
        with open(self.__correct) as correct_file:
            reader = Reader(correct_file)
            for i in reader.alleles():
                self.assertIsInstance(i, Record)
