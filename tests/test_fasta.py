#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import pyfaidx
import tempfile
import unittest
from bioformats.fasta import RandomSequence
from bioformats.fasta import Writer
from six import iteritems

path = os.path.dirname(__file__)
os.chdir(path)


class TestWriter(unittest.TestCase):
    def setUp(self):
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Check if sequences are written to a FASTA file properly.
        """
        output_file = self.__output_file

        # prepare the sequences
        fasta_patterns = ['AC', 'ACG', 'CCGT']
        sequences = {}
        for i, pattern in enumerate(fasta_patterns):
            sequences['chr{}'.format(i)] = pattern * 100

        with Writer(output_file) as output_fasta:
            for header, sequence in iteritems(sequences):
                output_fasta.write(header, sequence)

        # check the written file
        reader = pyfaidx.Fasta(output_file)
        for (header, sequence) in iteritems(sequences):
            self.assertEqual(sequence, reader[header][:].seq)

        os.unlink(output_file)
        os.unlink(output_file + '.fai')


class TestRandomSequence(unittest.TestCase):
    def test_get(self):
        """
        Check if a random sequence is generated.
        """
        seq_generator = RandomSequence(10)
        sequence = seq_generator.get()
        self.assertEqual(len(sequence), 10)
        self.assertIsInstance(sequence, str)
