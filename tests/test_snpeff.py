#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import unittest
from bioformats.snpeff import parse_snpeff_annotation, Record
from bioformats.snpeff import parse_hgvs_dna, HgvsRecord
from bioformats.exception import SnpEffError


class TestParseSnpEffAnnotation(unittest.TestCase):
    def setUp(self):
        # disable logging messages
        logging.disable(logging.ERROR)

    def test_parse_snpeff_annotation(self):
        line = 'ANN=C|missense_variant|MODERATE|OR4F5|' \
               'ENSG00000186092|transcript|ENST00000335137|' \
               'protein_coding|1/1|c.718G>C|p.Val240Leu|718/918|' \
               '718/918|240/305||'
        self.assertIsInstance(parse_snpeff_annotation(line), Record)

        with self.assertRaises(SnpEffError):
            parse_snpeff_annotation(line[4:])

        with self.assertRaises(SnpEffError):
            parse_snpeff_annotation(line[:-50])

    def test_parse_hgvs_dna(self):
        line = 'c.718G>C'
        self.assertIsInstance(parse_hgvs_dna(line), HgvsRecord)

        with self.assertRaises(SnpEffError):
            line = '718G>C'
            parse_hgvs_dna(line)

        with self.assertRaises(SnpEffError):
            line = 'c.718GC'
            parse_hgvs_dna(line)
            
        with self.assertRaises(SnpEffError):
            line = 'c.7G8G>C'
            parse_hgvs_dna(line)
