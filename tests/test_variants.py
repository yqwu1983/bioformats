#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import tempfile
import unittest
import bioformats.variants

path = os.path.dirname(__file__)
os.chdir(path)


class TestVcf2Genotypes(unittest.TestCase):
    def setUp(self):
        self.__input = os.path.join('data', 'variants', 'test.vcf.gz')
        self.__individuals = os.path.join('data', 'variants',
                                          'individuals.txt')
        self.__output = tempfile.NamedTemporaryFile().name

    def test_vcf2genotypes(self):
        with open(self.__input) as vcf_file:
            bioformats.variants.convert_vcf2genotypes(vcf_file,
                                                      self.__output)
        # now check with the list of individuals
        individuals = []
        with open(self.__individuals) as individuals_file:
            for line in individuals_file:
                individuals.append(line.rstrip())
        with open(self.__input) as vcf_file:
            bioformats.variants.convert_vcf2genotypes(vcf_file,
                                                      self.__output,
                                                      individuals)

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
