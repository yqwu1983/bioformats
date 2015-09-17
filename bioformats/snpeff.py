#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
from collections import namedtuple
from .exception import SnpEffError

logging.basicConfig()
logger = logging.getLogger(__name__)

snpeff_columns = ('allele', 'annotation', 'putative_impact',
                  'gene_name', 'gene_id', 'feature_type',
                  'feature_id', 'transcript_biotype', 'rank',
                  'total', 'hgvs_c', 'hgvs_p', 'cdna_pos',
                  'cdna_len', 'cds_pos', 'cds_len', 'protein_pos',
                  'protein_len', 'distance', 'extra')

Record = namedtuple('Record', snpeff_columns)

HgvsRecord = namedtuple('HgvsRecord', ('pos', 'ref', 'alt'))


def parse_hgvs_dna(x):
    """
    Parse a HGVS notation record.

    :param x: a DNA variant in the HGVS notation
    :type x: str
    :return: a parsed HGVS record
    :rtype: HgvsRecord
    """
    if not x.startswith('c.'):
        logger.error('incorrect HGVS DNA notation %s', x)
        raise SnpEffError
    x = x[2:]
    x_parts = x.split('>', 1)
    if len(x_parts) < 2:
        logger.error('incorrect HGVS DNA notation %s', x)
        raise SnpEffError
    pos = x_parts[0][:-1]
    try:
        pos = int(pos)
    except ValueError:
        logger.error('incorrect numeric value %s in HGVS notation',
                     pos)
        raise SnpEffError
    ref, alt = x_parts[0][-1], x_parts[1]

    return HgvsRecord(pos=pos, ref=ref, alt=alt)


def parse_snpeff_annotation(annotation):
    """
    Given a string representing an snpEff annotation, parse it.

    :param annotation: a string representing an snpEff annotation
    :type annotation: str
    :return: a parsed snpEff annotation record
    :rtype: Record
    """
    # check the annotation header
    if len(annotation) < 4 or annotation[:4] != 'ANN=':
        logger.error('incorrect snpEff annotation line header')
        raise SnpEffError
    annotation = annotation[4:]
    ann_parts = annotation.split('|', 15)
    if len(ann_parts) < 16:
        logger.error('incorrect number of snpEff annotation fields')
        raise SnpEffError
    record_fields = ann_parts[:8]
    # rank/total field
    record_fields += map(int, ann_parts[8].split('/'))
    # HGVS.c and HGVS.p fields
    record_fields += ann_parts[9:11]
    # cDNA position/length field
    record_fields += map(int, ann_parts[11].split('/'))
    # CDS position/length field
    record_fields += map(int, ann_parts[12].split('/'))
    # Protein position/length field
    record_fields += map(int, ann_parts[13].split('/'))
    # feature distance field
    record_fields += [ann_parts[14]]
    # errors, warnings of information messages field
    record_fields.append(ann_parts[15])

    return Record(*record_fields)
