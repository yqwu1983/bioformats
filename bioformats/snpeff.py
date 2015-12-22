#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import vcf
from collections import namedtuple
from . import bed
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


annotation_impacts = {
    'chromosome_number_variation': 'high',
    'exon_loss_variant': 'high',
    'frameshift_variant': 'high',
    'rare_amino_acid_variant': 'high',
    'splice_acceptor_variant': 'high',
    'splice_donor_variant': 'high',
    'start_lost': 'high',
    'stop_gained': 'high',
    'stop_lost': 'high',
    'transcript_ablation': 'high',
    '3_prime_UTR_truncation': 'moderate',
    '5_prime_UTR_truncation': 'moderate',
    'coding_sequence_variant': 'moderate',
    'disruptive_inframe_deletion': 'moderate',
    'disruptive_inframe_insertion': 'moderate',
    'inframe_deletion': 'moderate',
    'inframe_insertion': 'moderate',
    'missense_variant': 'moderate',
    'regulatory_region_ablation': 'moderate',
    'splice_region_variant': 'moderate',
    'TFBS_ablation': 'moderate',
    'sequence_feature': 'low',
    '5_prime_UTR_premature_start_codon_gain_variant': 'low',
    'initiator codon variant': 'low',
    'start_retained': 'low',
    'stop_retained variant': 'low',
    'synonymous variant': 'low',
    '3_prime_UTR_variant': 'modifier',
    '5_crime_UTR_variant': 'modifier',
    'conserved_intergenic_variant': 'modifier',
    'conserved_intron_variant': 'modifier',
    'downstream_gene_variant': 'modifier',
    'exon_variant': 'modifier',
    'feature_elongation': 'modifier',
    'feature_truncation': 'modifier',
    'gene_variant': 'modifier',
    'intergenic_region': 'modifier',
    'intragenic_variant': 'modifier',
    'intron_variant': 'modifier',
    'mature_miRNA_variant': 'modifier',
    'miRNA': 'modifier',
    'NMD_transcript_variant': 'modifier',
    'non_coding_transcript_exon_variant': 'modifier',
    'non_coding_transcript_variant': 'modifier',
    'regulatory_region_amplification': 'modifier',
    'regulatory_region_variant': 'modifier',
    'TF_binding_site_variant': 'modifier',
    'FFBS_amplification': 'modifier',
    'transcript_amplification': 'modifier',
    'transcript_variant': 'modifier',
    'upstream_gene_variant': 'modifier'
}

impact_colors = {
    'high': (255, 0, 0),
    'moderate': (0, 255, 255),
    'low': (0, 0, 255),
    'modifier': (255, 255, 255)
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
                        bed_color = impact_colors[annotation_impacts[
                            var_ann.annotation]]
                        bed_extra = list(var_ann)
                        if var_ann.hgvs_c is not None:
                            bed_extra[10] = var_ann.hgvs_c
                        else:
                            bed_extra[10] = 'NA'
                        if var_ann.hgvs_p is not None:
                            bed_extra[11] = 'p.{}{}>{}'.format(
                                var_ann.hgvs_p.pos, var_ann.hgvs_p.ref,
                                var_ann.hgvs_p.alt
                            )
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
