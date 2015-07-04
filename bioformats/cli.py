#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import pyfaidx
import re
from . import fasta


def randomfasta():
    """
    This function corresponds to a command-line tool which creates a
    FASTA file with random sequences of the specified number and
    length.
    """
    parser = argparse.ArgumentParser('Create a FASTA file with random '
                                     'nucleotide sequences')
    parser.add_argument('seq_length', type=int,
                        help='random sequence length')
    parser.add_argument('seq_num', type=int, help='random sequence '
                                                  'number')
    parser.add_argument('output', help='output filename')
    args = parser.parse_args()

    seq_generator = fasta.RandomSequence(args.seq_length)
    with fasta.Writer(args.output) as output_fasta:
        for i in range(args.seq_num):
            output_fasta.write('random_seq{}'.format(i+1),
                               seq_generator.get())


def fastagaps():
    """
    This functions corresponds to a command-line tool which detect
    gap regions in sequences of the specified FASTA file.
    :return:
    """
    parser = argparse.ArgumentParser(
        description='Get coordinates of gap regions in FASTA file '
                    'sequences')
    parser.add_argument('fasta_file', help='a FASTA file')
    parser.add_argument('bed_gaps', help='an output BED file of gap'
                                         'regions')
    args = parser.parse_args()

    gap_pattern = "[Nn-]+"
    with open(args.bed_gaps, 'w') as output:
        for seq in pyfaidx.Fasta(args.fasta_file):
            gap_regions = [(m.start(), m.end()) for m in re.finditer(
                gap_pattern, str(seq))]
            for start, end in gap_regions:
                output.write('{}\t{}\t{}\n'.format(
                    seq.name, start, end))
