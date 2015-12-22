#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import bioformats.bed
import logging
import unittest
import os
import tempfile
from bioformats.snpeff import parse_snpeff_ann, Record
from bioformats.snpeff import parse_hgvs_prot
from bioformats.snpeff import convert_snpeff2bed
from bioformats.snpeff import HgvsRecord
from bioformats.exception import SnpEffError

path = os.path.dirname(__file__)
os.chdir(path)


class TestParseSnpEffAnnotation(unittest.TestCase):
    def setUp(self):
        self.__vcf_file = os.path.join(
            'data', 'snpeff', 'snpeff.vcf'
        )
        self.__output = tempfile.NamedTemporaryFile().name
        # disable logging messages
        logging.disable(logging.ERROR)

    def test_parse_snpeff_annotation(self):
        line = 'C|missense_variant|MODERATE|OR4F5|' \
               'ENSG00000186092|transcript|ENST00000335137|' \
               'protein_coding|1/1|c.718G>C|p.Val240Leu|718/918|' \
               '718/918|240/305||'
        result = parse_snpeff_ann(line)
        self.assertIsInstance(result, Record)
        self.assertIsInstance(result.hgvs_p, HgvsRecord)

        line = 'C|missense_variant|MODERATE|OR4F5|' \
               'ENSG00000186092|transcript|ENST00000335137|' \
               'protein_coding||||||||'
        parse_snpeff_ann(line)

        with self.assertRaises(SnpEffError):
            parse_snpeff_ann(line[:-50])

    def test_parse_hgvs_prot(self):
        line = 'p.Val240Leu'
        result = parse_hgvs_prot(line)
        self.assertIsInstance(result, HgvsRecord)
        self.assertEqual(result.pos, 240)
        self.assertEqual(result.ref, 'Val')
        self.assertEqual(result.alt, 'Leu')

        with self.assertRaises(SnpEffError):
            line = 'Val240Leu'
            parse_hgvs_prot(line)

        with self.assertRaises(SnpEffError):
            line = 'p.ValLeu'
            parse_hgvs_prot(line)

        with self.assertRaises(SnpEffError):
            line = 'p.Val2E0Leu'
            parse_hgvs_prot(line)

    def test_convert_snpeff2bed(self):
        # convert to a BED file and try to read it
        convert_snpeff2bed(self.__vcf_file, self.__output)
        with open(self.__output) as bed_file:
            reader = bioformats.bed.Reader(bed_file)
            for record in reader.records():
                self.assertIsInstance(record, bioformats.bed.Record)
            self.assertEqual(reader.bed_columns, 9)
        # now check if a BED3 file is correctly generated
        convert_snpeff2bed(self.__vcf_file, self.__output, is_bed3=True)
        with open(self.__output) as bed_file:
            reader = bioformats.bed.Reader(bed_file)
            for record in reader.records():
                self.assertIsInstance(record, bioformats.bed.Record)

    def tearDown(self):
        if os.path.isfile(self.__output):
            os.unlink(self.__output)
