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

aa_code = {
    'Gly': 'G',
    'Pro': 'P',
    'Ala': 'A',
    'Val': 'V',
    'Leu': 'L',
    'Ile': 'I',
    'Met': 'M',
    'Cys': 'C',
    'Phe': 'F',
    'Tyr': 'Y',
    'Trp': 'W',
    'His': 'H',
    'Lys': 'K',
    'Arg': 'R',
    'Gln': 'Q',
    'Asn': 'N',
    'Glu': 'E',
    'Asp': 'D',
    'Ser': 'S',
    'Thr': 'T'
}


def parse_hgvs_dna(x):
    """
    Parse a HGVS DNA notation record.

    :param x: a DNA variant in the HGVS notation
    :type x: str
    :return: a parsed HGVS DNA record
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
        logger.error('incorrect position %s in HGVS notation', pos)
        raise SnpEffError
    ref, alt = x_parts[0][-1], x_parts[1]

    return HgvsRecord(pos=pos, ref=ref, alt=alt)


def parse_hgvs_prot(x):
    """
    Parse a HGVS protein notation record.

    :param x: a protein variant in the HGVS notation
    :type x: str
    :return: a parsed HGVS protein record
    :return: HgvsRecord
    """
    if not x.startswith('p.'):
        logger.error('incorrect HGVS protein notation %s', x)
        raise SnpEffError
    x = x[2:]
    if len(x) <= 6:
        logger.error('incorrect HGVS protein notation %s', x)
        raise SnpEffError
    ref = x[:3]
    alt = x[-3:]
    pos = x[3:-3]
    try:
        pos = int(pos)
    except ValueError:
        logger.error('incorrect position %s in HGVS notation', pos)
        raise SnpEffError

    return HgvsRecord(pos=pos, ref=ref, alt=alt)


def parse_snpeff_ann(annotation):
    """
    Given a string representing an snpEff annotation, parse it.

    :param annotation: a string representing an snpEff annotation
    :type annotation: str
    :return: a parsed snpEff annotation record
    :rtype: Record
    """
    ann_parts = annotation.split('|', 15)
    if len(ann_parts) < 16:
        logger.error('incorrect number of snpEff annotation fields')
        raise SnpEffError
    record_fields = ann_parts[:8]
    # rank/total field
    if ann_parts[8]:
        record_fields += map(int, ann_parts[8].split('/'))
    else:
        record_fields += [None, None]
    # HGVS.c and HGVS.p fields
    if ann_parts[9]:
        record_fields += [parse_hgvs_dna(ann_parts[9])]
    else:
        record_fields += [None]
    if ann_parts[10]:
        record_fields += [parse_hgvs_prot(ann_parts[10])]
    else:
        record_fields += [None]
    # cDNA position/length field
    if ann_parts[11]:
        record_fields += map(int, ann_parts[11].split('/'))
    else:
        record_fields += [None, None]
    # CDS position/length field
    if ann_parts[12]:
        record_fields += map(int, ann_parts[12].split('/'))
    else:
        record_fields += [None, None]
    # Protein position/length field
    if ann_parts[13]:
        record_fields += map(int, ann_parts[13].split('/'))
    else:
        record_fields += [None, None]
    # feature distance field
    record_fields += [ann_parts[14]]
    # errors, warnings of information messages field
    record_fields.append(ann_parts[15])

    return Record(*record_fields)
