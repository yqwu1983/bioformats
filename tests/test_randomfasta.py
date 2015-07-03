#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli


class TestRandomFasta(unittest.TestCase):
    def setUp(self):
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_randomfasta(self):
        """
        Check if a random FASTA file is created.
        """
        sys.argv = ['', '10', '10', self.__output_file]
        bioformats.cli.randomfasta()
        os.unlink(self.__output_file)
