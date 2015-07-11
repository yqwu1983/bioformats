#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import pyfaidx
import re
from . import fasta
from . import seqname


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


def renameseq():
    """
    This function corresponds to a command-line tool which renames
    sequences according to the specified renaming table.
    """
    parser = argparse.ArgumentParser(
        description='Change sequence names in a FASTA or plain text '
                    'file.'
    )

    # required parameters
    parser.add_argument('renaming_table',
                        help='a file containing a table of original '
                             'and new sequence names')
    parser.add_argument('input_file',
                        help='a file to change sequence names in')
    parser.add_argument('output_file',
                        help='an output file with renamed sequences')

    # optional arguments
    parser.add_argument('-f', '--fasta', action='store_true',
                        help='the input file is of the FASTA format')
    parser.add_argument('-c', '--column', type=int, default=1,
                        help='the number of the column that contains '
                             'sequence names to be changed staring '
                             'from 1')
    parser.add_argument('-r', '--revert', action='store_true',
                        help='perform reverse renaming, that is, '
                             'change original and new names in the '
                             'renaming table')
    parser.add_argument('--no_description', action='store_true',
                        help='remove descriptions from FASTA sequence '
                             'names')

    # parameters which specify the format of an input plain-text file
    parser.add_argument('--comment_char', default='#',
                        help='a character that designates comment '
                             'lines in the specified plain-text file')
    parser.add_argument('-s', '--separator', default='\t',
                        help='a symbol that separates columns in the '
                             'specified plain-text file')

    args = parser.parse_args()

    # according to the -f (--fasta) command-line option, choose the
    # appropriate renamer object
    if args.fasta:
        renamer = seqname.FastaSeqRenamer()
    else:
        renamer = seqname.TableSeqRenamer()

    renamer.read_renaming_dict(args.renaming_table)

    # for a user, the ecolumn numbers start from 1, for the program
    # from 0
    args.column -= 1

    if args.fasta:
        renamed_lines = renamer.renamed(args.input_file,
                                        reverse=args.revert)
    else:
        renamed_lines = renamer.renamed(args.input_file, args.column,
                                        sep=args.separator,
                                        comment_char=args.comment_char,
                                        reverse=args.revert)

    with open(args.output_file, 'w') as output:
        for line in renamed_lines:
            # if we are processing a FASTA file and the option to remove
            # descriptions from FASTA sequence headers is specified,
            # then remove the description
            if args.fasta and args.no_description:
                line = line.split()[0]
            output.write(line)
            if not args.fasta:
                output.write('\n')
