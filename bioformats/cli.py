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
                    'tabular file.'
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
                                        reverse=args.revert,
                                        no_desc=args.no_description)
    else:
        renamed_lines = renamer.renamed(args.input_file, args.column,
                                        sep=args.separator,
                                        comment_char=args.comment_char,
                                        reverse=args.revert)

    with open(args.output_file, 'w') as output:
        for line in renamed_lines:
            output.write(line)


def ncbirenameseq():
    """
    This function corresponds to a command-line tool that changes
    NCBI-style identificators of sequences in a tabular or FASTA file.
    """
    parser = argparse.ArgumentParser(
        description='Change NCBI-style sequence names in a FASTA file'
                    'or plain text tabular file',
        epilog='Format values: refseq_full, genbank_full, refseq_gi, '
               'genbank_gi, refseq, genbank, chr_refseq, chr_genbank '
               'and ucsc (only as the output format).'
    )

    # required parameters
    parser.add_argument('input_file',
                        help='a file to change sequence names in')
    parser.add_argument('input_fmt',
                        help='a format of sequence names in input')
    parser.add_argument('output_file',
                        help='an output file for renamed sequences')
    parser.add_argument('output_fmt',
                        help='a format of sequence names in output')

    # the input file format
    parser.add_argument('-f', '--fasta', action='store_true',
                        help='the input file is of the FASTA format')
    parser.add_argument('-c', '--column', type=int, default=1,
                        help='the number of the column that contains '
                             'sequence names to be changed (1 by '
                             'default)')

    # the format of an input plain text file
    parser.add_argument('--comment_char', default='#',
                        help='a character designating comment lines '
                             'in the specified plain text file')
    parser.add_argument('-s', '--separator', default='\t',
                        help='a symbol separating columns in the '
                             'specified plain text file')
    
    # NCBI accession number files
    parser.add_argument('--chr', 
                        help='a name of a file containing NCBI '
                             'chromosome accession numbers')
    parser.add_argument('--unloc', 
                        help='a name of a file containing NCBI '
                             'accession numbers of unlocalized '
                             'fragments')
    parser.add_argument('--unpl',
                        help='a name of a file containing NCBI '
                             'accession numbers of unplaced fragments')
    
    # prefixes for sequence names
    parser.add_argument('--prefix', default='',
                        help='a prefix to be added to sequence names')
    parser.add_argument('--prefix_chr', default='',
                        help='a prefix to be added to chromosome '
                             'names')
    parser.add_argument('--prefix_unloc', default='',
                        help='a prefix to be added to unlocalized '
                             'fragment names')
    parser.add_argument('--prefix_unpl', default='',
                        help='a prefix to be added to unplaced '
                             'fragment names')
    
    # suffixes for sequence names
    parser.add_argument('--suffix', default='',
                        help='a suffix to be added to sequence names')
    parser.add_argument('--suffix_chr', default='',
                        help='a suffix to be added to chromosome '
                             'names')
    parser.add_argument('--suffix_unloc', default='',
                        help='a suffix to be added to unlocalized '
                             'fragment names')
    parser.add_argument('--suffix_unpl', default='',
                        help='a suffix to be added to unplaced '
                             'fragment names')
    
    # auxiliary options for sequence names
    parser.add_argument('-r', '--revert', action='store_true',
                        help='perform reverse renaming, that is, '
                             'change original and new names in the '
                             'renaming table')
    parser.add_argument('--no_version', action='store_true',
                        help='remove a sequence version from an '
                             'accession number')
    parser.add_argument('--no_description', action='store_true',
                        help='remove descriptions from FASTA sequence '
                             'headers')
    parser.add_argument('--output_table',
                        help='write the renaming table to the '
                             'specified file')
    
    args = parser.parse_args()
    
    # choose the renaming object according to the specified 
    # command-line option
    if args.fasta:
        renamer = seqname.NcbiFastaSeqRenamer()
    else:
        renamer = seqname.NcbiTableSeqRenamer()

    if args.output_fmt == 'ucsc':
        chr_output_fmt = 'chr'
        unloc_output_fmt = unpl_output_fmt = 'chr_genbank'
        args.prefix = 'chr'
        args.prefix_chr = args.prefix_unloc = args.prefix_unpl = ''
        args.suffix = args.suffix_chr = args.suffix_unpl = ''
        args.suffix_unloc = '_random'
        args.no_version = True
    else:
        chr_output_fmt = unloc_output_fmt = unpl_output_fmt = \
            args.output_fmt

    # read accession numbers for chromosomes, unlocalized and 
    # unplaced fragments
    if args.chr:
        renamer.read_ncbi_acc_num(
            args.chr, args.input_fmt, chr_output_fmt,
            prefix=args.prefix + args.prefix_chr,
            suffix=args.suffix + args.suffix_chr,
            remove_seq_version=args.no_version
        )
    if args.unloc:
        renamer.read_ncbi_acc_num(
            args.unloc, args.input_fmt, unloc_output_fmt,
            prefix=args.prefix + args.prefix_unloc,
            suffix=args.suffix + args.suffix_unloc,
            remove_seq_version=args.no_version
        )
    if args.unpl:
        renamer.read_ncbi_acc_num(
            args.unpl, args.input_fmt, unpl_output_fmt,
            prefix=args.prefix + args.prefix_unpl,
            suffix=args.suffix + args.suffix_unpl,
            remove_seq_version=args.no_version
        )

    # for a user, the column numbers start from 1, but for the
    # program they start from 0
    args.column -= 1

    if args.fasta:
        renamed_lines = renamer.renamed(args.input_file,
                                        no_desc=args.no_description,
                                        reverse=args.revert)
    else:
        renamed_lines = renamer.renamed(args.input_file, args.column,
                                        sep=args.separator,
                                        comment_char=args.comment_char,
                                        reverse=args.revert)

    with open(args.output_file, 'w') as output_file:
        for line in renamed_lines:
            output_file.write(line)

    # write the renaming table to the specified file, if required
    if args.output_table:
        renamer.write_renaming_dict(args.output_table)
