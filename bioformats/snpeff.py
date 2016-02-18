#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import vcf
from collections import namedtuple
from collections import OrderedDict
from . import bed
from .variants import allele_pair_iterator
from .exception import SnpEffError
from future.utils import iteritems

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

impact_colors = {
    'HIGH': (255, 0, 0),
    'MODERATE': (0, 255, 255),
    'LOW': (0, 0, 255),
    'MODIFIER': (255, 255, 255)
}


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
        record_fields += [ann_parts[9]]
    else:
        record_fields += [None]
    if ann_parts[10]:
        if ann_parts[1] == 'missense_variant':
            record_fields += [parse_hgvs_prot(ann_parts[10])]
        else:
            record_fields += [ann_parts[10]]
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


def convert_snpeff2bed(vcf_file, bed_file, is_bed3=False):
    """
    Convert a specified snpEff-annotated VCF file to the BED format
    considering its attributes.

    :param vcf_file: a name of an snpEff-annotated VCF file
    :param bed_file: a name of the output BED file
    :param is_bed3: convert to the BED3 format
    :type vcf_file: str
    :type bed_file: str
    :type is_bed3: bool
    """
    total_processed = 0
    with open(vcf_file) as input_file:
        vcf_reader = vcf.Reader(input_file)
        with bed.Writer(bed_file) as bed_writer:
            for variant in vcf_reader:
                if 'ANN' in variant.INFO:
                    for var_ann in variant.INFO['ANN']:
                        var_ann = parse_snpeff_ann(var_ann)
                        bed_size = len(var_ann.allele) - \
                            len(variant.REF) + 1
                        bed_start = variant.POS - 1
                        bed_end = bed_start + bed_size
                        bed_color = impact_colors[
                            var_ann.putative_impact]
                        bed_extra = list(var_ann)
                        if var_ann.hgvs_c is not None:
                            bed_extra[10] = var_ann.hgvs_c
                        else:
                            bed_extra[10] = 'NA'
                        if var_ann.hgvs_p is not None:
                            if var_ann.annotation == \
                                    'missense_variant':
                                bed_extra[11] = 'p.{}{}>{}'.format(
                                    var_ann.hgvs_p.pos,
                                    var_ann.hgvs_p.ref,
                                    var_ann.hgvs_p.alt
                                )
                            else:
                                bed_extra[11] = var_ann.hgvs_p
                        else:
                            bed_extra[11] = 'NA'
                        if is_bed3:
                            bed_record = bed.Record(
                                seq=variant.CHROM,
                                start=bed_start,
                                end=bed_end,
                                name=None,
                                score=None,
                                strand=None,
                                thick_start=None,
                                thick_end=None,
                                color=None,
                                block_num=None,
                                block_sizes=None,
                                block_starts=None,
                                extra=bed_extra
                            )
                        else:
                            bed_record = bed.Record(
                                seq=variant.CHROM,
                                start=bed_start,
                                end=bed_end,
                                name=var_ann.annotation,
                                score=1000,
                                strand='+',
                                thick_start=bed_start,
                                thick_end=bed_end,
                                color=','.join(map(str, bed_color)),
                                block_num=None,
                                block_sizes=None,
                                block_starts=None,
                                extra=bed_extra
                            )
                        bed_writer.write(bed_record)
                        total_processed += 1
    logger.info('%d variant effects processed', total_processed)


def gene_feature_id_iterator(record):
    """
    Given a record from an snpEff-annotated VCF file, iterate over
    pairs of gene and feature IDs of its effects.

    :param record: a record from an snpEff-annotated VCF file
    :type record: vcf.model._Record
    :return: an iterator to iterate over pairs of gene and feature
        IDs of variant effects
    """
    if 'ANN' in record.INFO:
        gene_feature_pairs = set()
        for effect in record.INFO['ANN']:
            ann = parse_snpeff_ann(effect)
            gene_feature_pairs.add((ann.gene_id, ann.feature_id))
        for i in sorted(list(gene_feature_pairs)):
            yield i
    else:
        yield


def sample_effect_iterator(record):
    """
    Given a record from an snpEff-annotated VCF file, iterate over
    effects of the related genotypes.

    :param record: a record from an snpEff-annotated VCF file
    :type record: vcf.model._Record
    :return: an iterator to iterate over pairs of gene and feature
        IDs of variant effects
    """
    if 'ANN' in record.INFO:
        # first, form the dictionary of variant effects which keys
        # are tuples of three elements: an alternative allele,
        # a gene ID and a transcript ID
        effects = {}
        for ann in record.INFO['ANN']:
            ann = parse_snpeff_ann(ann)
            key = (ann.allele, ann.gene_id, ann.feature_id)
            if key not in effects:
                effects[key] = ann.annotation
            elif effects[key] == ann.annotation:
                logger.warning('duplicate effect record: (%s,%s,'
                               '%s) - %s',
                               ann.allele,
                               ann.gene_id,
                               ann.feature_id,
                               effects[key])
            else:
                logger.error('conflicting effect records: (%s,%s,'
                             '%s) - %s vs %s !',
                             ann.allele,
                             ann.gene_id,
                             ann.feature_id,
                             effects[key],
                             ann.effect)
                raise SnpEffError
        # next, iterate over genotypes and return tuples of two
        # elements: an original genotype and its effects
        gene_feature_ids = list(gene_feature_id_iterator(record))
        genotypes_dict = OrderedDict()
        for i in allele_pair_iterator(record):
            genotypes_dict[i] = OrderedDict()
            for g in gene_feature_ids:
                genotypes_dict[i][g] = dict()
        for g in record.samples:
            alleles = tuple(sorted(g.gt_bases.split(g.gt_phase_char())))
            for e in effects:
                genotypes_dict[alleles][(e[1], e[2])] = (
                    g.sample,
                    effects[(alleles[0], e[1], e[2])]
                    if alleles[0] != record.REF else '-',
                    effects[(alleles[1], e[1], e[2])]
                    if alleles[1] != record.REF else '-'
                )
        # we have formed the dictionary of effects in the required
        # order; now release its values in that order
        for i, j in iteritems(genotypes_dict):
            for k, l in iteritems(j):
                yield i + k + l
    else:
        yield
