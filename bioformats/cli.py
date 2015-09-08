#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import pyfaidx
import re
import sys
from . import autosql
from . import fasta
from . import seqname
from . import bed
from . import exception
from argparse import RawTextHelpFormatter


def bioformats():
    """
    The main function that is run if bioformats was launched. It
    processes command-line arguments to be passed to subroutines.
    """
    parser = argparse.ArgumentParser(
        description='A collection of tools to process bioinformatic '
                    'data.'
    )

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 0.1.5')

    subparsers = parser.add_subparsers(dest='command')

    randomfasta_parser(subparsers)
    fastagaps_parser(subparsers)
    renameseq_parser(subparsers)
    ncbirenameseq_parser(subparsers)
    fastareorder_parser(subparsers)
    bedcolumns_parser(subparsers)
    bedautosql_parser(subparsers)

    args = parser.parse_args()

    launchers = dict([
        ('randomfasta', randomfasta_launcher),
        ('fastagaps', fastagaps_launcher),
        ('renameseq', renameseq_launcher),
        ('ncbirenameseq', ncbirenameseq_launcher),
        ('fastareorder', fastareorder_launcher),
        ('bedcolumns', bedcolumns_launcher),
        ('bedautosql', bedautosql_launcher)
    ])

    launchers[args.command](args)


def randomfasta_parser(subparsers):
    """
    Parser for the randomfasta tool.
    """
    parser = subparsers.add_parser(
        'randomfasta',
        help='create a random FASTA file',
        description='Create a FASTA file with random nucleotide '
                    'sequences.'
    )
    parser.add_argument('seq_length', type=int,
                        help='random sequence length')
    parser.add_argument('seq_num', type=int, help='random sequence '
                                                  'number')
    parser.add_argument('output', help='output filename')


def randomfasta_launcher(args):
    """
    Launcher for the randomfasta tool.
    """
    seq_generator = fasta.RandomSequence(args.seq_length)
    with fasta.Writer(args.output) as output_fasta:
        for i in range(args.seq_num):
            output_fasta.write('random_seq{}'.format(i+1),
                               seq_generator.get())


def fastagaps_parser(subparsers):
    """
    Parser for the fastagaps tool.
    """
    parser = subparsers.add_parser(
        'fastagaps',
        help='get gaps from a FASTA file',
        description='Get coordinates of gap regions in sequences of '
                    'a FASTA file.'
    )
    parser.add_argument('fasta_file', help='a FASTA file')
    parser.add_argument('bed_gaps', help='an output BED file of gap'
                                         'regions')


def fastagaps_launcher(args):
    """
    Launcher for the fastagaps tool.
    """
    gap_pattern = "[Nn-]+"
    with open(args.bed_gaps, 'w') as output:
        for seq in pyfaidx.Fasta(args.fasta_file):
            gap_regions = [(m.start(), m.end()) for m in re.finditer(
                gap_pattern, str(seq))]
            for start, end in gap_regions:
                output.write('{}\t{}\t{}\n'.format(
                    seq.name, start, end))


def renameseq_parser(subparsers):
    """
    Parser for the renameseq tool.
    """
    parser = subparsers.add_parser(
        'renameseq',
        help='rename sequences in a FASTA or tabular file',
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


def renameseq_launcher(args):
    """
    Launcher for the renameseq tool.
    """
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


def ncbirenameseq_parser(subparsers):
    """
    Parser for the ncbirenameseq tool.
    """
    parser = subparsers.add_parser(
        'ncbirenameseq',
        help='rename NCBI-named sequences in a FASTA if tabular file',
        description='Change NCBI-style sequence names in a FASTA file'
                    'or plain text tabular file',
        epilog='Format values: \n'
               '\trefseq_full:\tgi|568815597|ref|NC_000001.11|\n'
               '\tgenbank_full:\tgi|568336023|gb|CM000663.2|\n'
               '\trefseq_gi:\t568815597\n'
               '\tgenbank_gi:\t568336023\n'
               '\trefseq: \tNC_000001.11\n'
               '\tgenbank:\tCM000663.2\n'
               '\tchr_refseq:\t1_NC_000001.11\n'
               '\tchr_genbank:\t1_CM000663.2',
        formatter_class=RawTextHelpFormatter
    )

    format_values = ['refseq_full', 'genbank_full', 'refseq_gi',
                     'genbank_gi', 'refseq', 'genbank', 'chr_refseq',
                     'chr_genbank']

    # required parameters
    parser.add_argument('input_file',
                        help='a file to change sequence names in')
    parser.add_argument('input_fmt',
                        choices=format_values,
                        help='a format of sequence names in input')
    parser.add_argument('output_file',
                        help='an output file for renamed sequences')
    parser.add_argument('output_fmt',
                        choices=format_values + ['ucsc'],
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


def ncbirenameseq_launcher(args):
    """
    Launcher for the ncbirenameseq tool.
    """
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


def fastareorder_parser(subparsers):
    """
    Parser for the fastareorder tool.
    """
    parser = subparsers.add_parser(
        'fastareorder',
        help='reorder sequences in a FASTA file',
        description='Reorder sequences in a FASTA file.'
    )

    parser.add_argument('fasta', help='a FASTA file of sequences to '
                                      'reorder')
    parser.add_argument('order_file', help='a file with the sequence '
                                           'order')
    parser.add_argument('output', help='an output FASTA file of '
                                       'reordered sequences')

    parser.add_argument('-i', '--ignore_missing', action='store_true',
                        help='ignore sequences in the specified order '
                             'file that are missing in the input '
                             'FASTA file')


def fastareorder_launcher(args):
    """
    Launcher for the fastareorder tool.
    """
    reorderer = fasta.Reorder(args.order_file)
    reorderer.write(args.fasta, args.output,
                    ignore_missing=args.ignore_missing)


def bedcolumns_parser(subparsers):
    """
    Parser for the bedcolumns tool.
    """
    parser = subparsers.add_parser(
        'bedcolumns',
        help='determine the numbers of BED and extra columns',
        description='Determine the number of BED columns and the '
                    'number of extra columns in a BED file'
    )

    parser.add_argument('bed_file', help='a BED file')


def bedcolumns_launcher(args):
    """
    Launcher for the bedcolumns tool.
    """
    with open(args.bed_file) as bed_file:
        reader = bed.Reader(bed_file)
        try:
            for _ in reader.records():
                pass
            if reader.aux_columns > 0:
                print('{}+{}'.format(reader.bed_columns,
                                     reader.aux_columns))
            else:
                print('{}'.format(reader.bed_columns))
        except exception.BedError:
            sys.stderr.write('Incorrect BED file.\n')


def bedautosql_parser(subparsers):
    """
    Parser for the bedautosql tool.
    """
    parser = subparsers.add_parser(
        'bedautosql',
        help='get an autoSql table for a BED file',
        description='Get an autoSql table structure for the specified '
                    'BED file'
    )

    parser.add_argument('bed_file', help='a BED file')
    parser.add_argument('output_file', help='an output file')

    # optional arguments
    parser.add_argument('-n', '--name', default='Table',
                        help='a table name')
    parser.add_argument('-d', '--description', default='Description',
                        help='a table description')


def bedautosql_launcher(args):
    """
    Launcher for the bedautosql tool.
    """
    with open(args.bed_file) as bed_file:
        reader = bed.Reader(bed_file)
        try:
            table = bed.get_autosql_table(reader, args.name,
                                          args.description)
            with autosql.Writer(args.output_file, table.name,
                                table.desc) as writer:
                for i in table.entries:
                    writer.write(i)
        except exception.BedError:
                sys.stderr.write('Incorrect BED file.\n')
