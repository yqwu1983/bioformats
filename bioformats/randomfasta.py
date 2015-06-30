#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
from . import fasta


def main():
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
