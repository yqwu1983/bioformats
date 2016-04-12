#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import os
import unittest
from bioformats.agp import Reader


path = os.path.dirname(__file__)
os.chdir(path)


class TestReader(unittest.TestCase):
    def setUp(self):
        self.correct_file = os.path.join(
            'data', 'agp', 'correct.agp'
        )

        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_agp(self):
        """
        Check if the parser reads an AGP file in the correct way.
        """
        with open(self.correct_file) as input_file:
            parser = Reader(input_file)
            for _ in parser.records():
                pass
