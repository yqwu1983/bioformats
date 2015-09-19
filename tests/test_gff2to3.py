#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import sys
import tempfile
import unittest
import bioformats.cli

path = os.path.dirname(__file__)
os.chdir(path)


class TestGff2to3(unittest.TestCase):
    def setUp(self):
        self.__input_file = os.path.join('data', 'gff3',
                                         'example.gff2')
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_gff2to3(self):
        sys.argv = ['', 'gff2to3', self.__input_file,
                    self.__output_file]
        bioformats.cli.bioformats()

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)
