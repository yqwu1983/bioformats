#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import pyfaidx
import tempfile
import unittest
from bioformats.bed import Reader
from bioformats.exception import BioformatsError
from bioformats.fasta import RandomSequence
from bioformats.fasta import Reorder
from bioformats.fasta import Writer
from bioformats.fasta import FlankNFilter
from future.utils import iteritems

try:
    import itertools.izip as zip
except ImportError:
    pass

path = os.path.dirname(__file__)
os.chdir(path)


class TestWriter(unittest.TestCase):
    def setUp(self):
        self.__output_file = tempfile.NamedTemporaryFile().name

    def test_write(self):
        """
        Check if sequences are written to a FASTA file properly.
        """
        # prepare the sequences
        fasta_patterns = ['AC', 'ACG', 'CCGT']
        sequences = {}
        for i, pattern in enumerate(fasta_patterns):
            sequences['chr{}'.format(i)] = pattern * 100

        with Writer(self.__output_file) as output_fasta:
            for header, sequence in iteritems(sequences):
                output_fasta.write(header, sequence)

        # check the written file
        reader = pyfaidx.Fasta(self.__output_file)
        for (header, sequence) in iteritems(sequences):
            self.assertEqual(sequence, reader[header][:].seq)

    def tearDown(self):
        for i in (self.__output_file, self.__output_file + '.fai'):
            if os.path.isfile(i):
                os.unlink(i)


class TestRandomSequence(unittest.TestCase):
    def test_get(self):
        """
        Check if a random sequence is generated.
        """
        seq_generator = RandomSequence(10)
        sequence = seq_generator.get()
        self.assertEqual(len(sequence), 10)
        self.assertIsInstance(sequence, str)


class TestReorder(unittest.TestCase):
    def setUp(self):
        self.__input = os.path.join(
            'data', 'fasta', 'test.fa'
        )
        self.__order = os.path.join(
            'data', 'fasta', 'order.txt'
        )
        self.__output = tempfile.NamedTemporaryFile().name

    def test_reorder(self):
        """
        Check if input sequences are properly reordered.
        """
        test = Reorder(self.__order)
        test.write(self.__input, self.__output, ignore_missing=True)
        # check if sequences are in the specified order
        input_fasta = pyfaidx.Fasta(self.__input)
        output_fasta = pyfaidx.Fasta(self.__output)
        present_seq = [x for x in test.order if x in input_fasta.keys()]
        self.assertEqual(present_seq, list(output_fasta.keys()))

        with self.assertRaises(BioformatsError):
            test.write(self.__input, self.__output,
                       ignore_missing=False)

    def tearDown(self):
        for i in (self.__input + '.fai',
                  self.__output,
                  self.__output + '.fai'):
            if os.path.isfile(i):
                os.unlink(i)


class TestFlankNFilter(unittest.TestCase):
    def setUp(self):
        self.test_fa = os.path.join(
            'data', 'fasta', 'flanknfilter_test.fa'
        )
        self.test_bed = os.path.join(
            'data', 'fasta', 'flanknfilter_test.bed'
        )
        self.output_results = dict([
            ((2, False), 'fnf_test_2f.bed'),
            ((3, False), 'fnf_test_3f.bed'),
            ((4, False), 'fnf_test_4f.bed'),
            ((4, True), 'fnf_test_4t.bed')
        ])
        self.output = tempfile.NamedTemporaryFile().name

    def test_bed(self):
        for i, j in iteritems(self.output_results):
            feature_filter = FlankNFilter(self.test_fa, flank_len=i[0])
            test_output = os.path.join('data', 'fasta', j)
            feature_filter.filter_bed(self.test_bed, self.output,
                                      i[1])
            with open(test_output) as correct_file:
                with open(self.output) as produced_file:
                    for k, l in zip(
                            Reader(correct_file).records(),
                            Reader(produced_file).records()):
                        self.assertEqual(k, l)

    def tearDown(self):
        if os.path.isfile(self.output):
            os.unlink(self.output)
