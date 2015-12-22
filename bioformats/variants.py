#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import operator
import vcf

logging.basicConfig()
logger = logging.getLogger(__name__)


def get_unique_genotypes(record, individuals=None):
    """
    Given a record from a VCF file, return a set of all unique
    genotypes specific to it.

    :param record: a record from a VCF file
    :param individuals: a list of names of individuals which
        genotypes are to be considered
    :type record: vcf.Record
    :type individuals: list
    :return: a set of its unique genotypes
    :rtype: set
    """
    if individuals is None:
        # consider all individuals from the file
        individuals = map(operator.attrgetter('sample'), record.samples)
    genotypes = set()
    for i in record.samples:
        if i.sample in individuals:
            new_genotype = sorted(i.gt_bases.split(i.gt_phase_char()))
            if len(new_genotype) < 2:
                logger.warning('the genotype %s at %s:%d contains '
                               'of the sample %s only a single '
                               'allele',
                               i.gt_bases,
                               record.CHROM,
                               record.POS,
                               i.sample)
                new_genotype.append('-')
            genotypes.add(tuple(new_genotype))
    return genotypes


def get_all_genotypes(record, individuals=None):
    """
    Given a record from a VCF file, return a list of all its genotypes.

    :param record: a record from a VCF file
    :param individuals: a list of names of individuals which
        genotypes are to be considered
    :type record: vcf.Record
    :type individuals: list
    :return: a list of all genotypes
    :rtype: list
    """
    if individuals is None:
        # consider all individuals from the file
        individuals = map(operator.attrgetter('sample'), record.samples)
    genotypes = []
    for i in record.samples:
        if i.sample in individuals:
            new_genotype = sorted(i.gt_bases.split(i.gt_phase_char()))
            if len(new_genotype) < 2:
                logger.warning('the genotype %s at %s:%d contains '
                               'of the sample %s only a single '
                               'allele',
                               i.gt_bases,
                               record.CHROM,
                               record.POS,
                               i.sample)
                new_genotype.append('-')
            genotypes.append(tuple(new_genotype))
    return genotypes


def convert_vcf2genotypes(vcf_file, output_filename, individuals=None):
    """
    Given a VCF file, prepare the genotype statistics file from it.

    :param vcf_file: a VCF file handle
    :param individuals: a list of individuals which genotypes are to
        be counted
    :param output_filename: an output file name
    :type output_filename: str
    :type individuals: list
    """
    template = '\t'.join(['{}'] * 9) + '\n'
    with open(output_filename, 'w') as output:
        for record in vcf.Reader(vcf_file):
            genotypes = get_all_genotypes(record, individuals)
            for i in get_unique_genotypes(record, individuals):
                hom_1 = hom_2 = het = 0
                cur_alleles = set(i)
                for j in genotypes:
                    if not cur_alleles.issuperset(set(j)):
                        continue
                    if j[1] == '-':
                        # we have only one allele for the variant
                        # it cannot be a heterozygote
                        if j[0] == i[0]:
                            hom_1 += 1
                    else:
                        # consider this genotype because it contains
                        # at least allele of the current genotype
                        if j[0] == j[1]:
                            if i[0] == j[0]:
                                hom_1 += 1
                            else:
                                hom_2 += 1
                        else:
                            het += 1
                if i[0] == record.REF:
                    ref_allele = 1
                elif i[1] == record.REF:
                    ref_allele = 2
                else:
                    ref_allele = 0
                output.write(template.format(
                    record.CHROM,
                    record.POS - 1,
                    record.POS,
                    i[0],
                    i[1],
                    hom_1,
                    het,
                    hom_2,
                    ref_allele
                ))
