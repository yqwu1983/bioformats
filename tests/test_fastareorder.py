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


class TestFastaReorder(unittest.TestCase):
    def setUp(self):
        self.__fasta = os.path.join(
            'data', 'fasta', 'test.fa'
        )
        self.__order = os.path.join(
            'data', 'fasta', 'order.txt'
        )
        self.__output = tempfile.NamedTemporaryFile().name

    def test_fastareorder(self):
        """
        Test if the tool runs correctly.
        """
        sys.argv = ['', 'fastareorder', '-i', self.__fasta,
                    self.__order, self.__output]
        bioformats.cli.bioformats()

    def tearDown(self):
        for i in (self.__fasta + '.fai',
                  self.__output,
                  self.__output + '.fai'):
            if os.path.isfile(i):
                os.unlink(i)
