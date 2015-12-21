#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import operator
import vcf


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
            genotypes.add(tuple(sorted(i.gt_bases.split(
                i.gt_phase_char()))))
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
            genotypes.append(tuple(sorted(i.gt_bases.split(
                i.gt_phase_char()))))
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
                if i[0] == i[1]:
                    # skip homozygous genotypes
                    continue
                hom_1 = hom_2 = het = 0
                for j in genotypes:
                    if i[0] in j or i[1] in j:
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
