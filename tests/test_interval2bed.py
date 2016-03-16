#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestInterval2Bed(unittest.TestCase):
    def setUp(self):
        self.interval = os.path.join('data', 'interval', 'correct.txt')
        self.output = tempfile.NamedTemporaryFile().name

    def test_interval2bed(self):
        sys.argv = ['', 'interval2bed', self.interval, self.output]
        bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.exists(self.output):
            os.unlink(self.output)
