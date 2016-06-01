#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import tempfile
import unittest
import bioformats.variants as variants
import bioformats.bed
import vcf

path = os.path.dirname(__file__)
os.chdir(path)


class TestVariantRoutines(unittest.TestCase):
    def setUp(self):
        self.__input = os.path.join('data', 'variants', 'test.vcf')
        self.__multiallele_input = os.path.join(
            'data', 'variants', 'multiallele_indels.vcf')
        self.__individuals = os.path.join('data', 'variants',
                                          'individuals.txt')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_vcf2genotypes(self):
        with open(self.__input) as vcf_file:
            variants.convert_vcf2genotypes(vcf_file,
                                                      self.__output)
        # now check with the list of individuals
        individuals = []
        with open(self.__individuals) as individuals_file:
            for line in individuals_file:
                individuals.append(line.rstrip())
        with open(self.__input) as vcf_file:
            variants.convert_vcf2genotypes(vcf_file,
                                                      self.__output,
                                                      individuals)

    def test_allele_pair_iterator(self):
        for filename in (self.__input, self.__multiallele_input):
            with open(filename) as vcf_file:
                reader = vcf.Reader(vcf_file)
                for v in reader:
                    for _ in variants.allele_pair_iterator(v):
                        pass

    def test_vcf2bed(self):
        variants.vcf2bed(self.__input, self.__output)
        # verify the obtained BED file
        with open(self.__output) as bed_file:
            bed_reader = bioformats.bed.Reader(bed_file)
            for record in bed_reader.records():
                self.assertIsInstance(record, bioformats.bed.Record)

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
