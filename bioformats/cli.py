#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import pyfaidx
import re
import sys
import vcf
from . import autosql
from . import fasta
from . import seqname
from . import bed
from . import gff3
from . import exception
from . import repeatmasker
from . import snpeff
from . import variants
from . import interval
from . import __version__


def bioformats():
    """
    The main function that is run if bioformats was launched. It
    processes command-line arguments to be passed to subroutines.
    """
    parser = argparse.ArgumentParser(
        description='A collection of tools to process bioinformatic '
                    'data.',
        usage='%(prog)s [command] [-h] [-v]\nPlease specify the '
              'command or use -h to view the help message.'
    )

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(dest='command')

    subparser_routines = {
        'randomfasta': randomfasta_parser,
        'fastagaps': fastagaps_parser,
        'renameseq': renameseq_parser,
        'ncbirenameseq': ncbirenameseq_parser,
        'fastareorder': fastareorder_parser,
        'bedcolumns': bedcolumns_parser,
        'bedautosql': bedautosql_parser,
        'rmout2bed': rmout2bed_parser,
        'gfftagstat': gfftagstat_parser,
        'gff2to3': gff2to3_parser,
        'snpeff2pph': snpeff2pph_parser,
        'gff2bed': gff2bed_parser,
        'snpeff2bed': snpeff2bed_parser,
        'vcfgeno2bed': vcfgeno2bed_parser,
        'vcfeffect2bed': vcfeffect2bed_parser,
        'flanknfilter': flanknfilter_parser,
        'interval2bed': interval2bed_parser,
        'vcf2bed': vcf2bed_parser
    }

    for i in sorted(subparser_routines):
        subparser_routines[i](subparsers)

    args = parser.parse_args()

    launchers = dict([
        ('randomfasta', randomfasta_launcher),
        ('fastagaps', fastagaps_launcher),
        ('renameseq', renameseq_launcher),
        ('ncbirenameseq', ncbirenameseq_launcher),
        ('fastareorder', fastareorder_launcher),
        ('bedcolumns', bedcolumns_launcher),
        ('bedautosql', bedautosql_launcher),
        ('rmout2bed', rmout2bed_launcher),
        ('gfftagstat', gfftagstat_launcher),
        ('gff2to3', gff2to3_launcher),
        ('snpeff2pph', snpeff2pph_launcher),
        ('gff2bed', gff2bed_launcher),
        ('snpeff2bed', snpeff2bed_launcher),
        ('vcfgeno2bed', vcfgeno2bed_launcher),
        ('vcfeffect2bed', vcfeffect2bed_launcher),
        ('flanknfilter', flanknfilter_launcher),
        ('interval2bed', interval2bed_launcher),
        ('vcf2bed', vcf2bed_launcher)
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
        formatter_class=argparse.RawTextHelpFormatter
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
                    'BED file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('bed_file', help='a BED file')
    parser.add_argument('output_file', help='an output file')

    # optional arguments
    parser.add_argument('-n', '--name', default='Table',
                        help='a table name')
    parser.add_argument('-d', '--description', default='Description',
                        help='a table description')
    parser.add_argument('-l', '--lines', type=int, default=100,
                        help='the number of lines to analyze'
                             'from the input file')


def bedautosql_launcher(args):
    """
    Launcher for the bedautosql tool.
    """
    with open(args.bed_file) as bed_file:
        reader = bed.Reader(bed_file)
        try:
            table = bed.get_autosql_table(reader, args.name,
                                          args.description,
                                          args.lines)
            with autosql.Writer(args.output_file, table.name,
                                table.desc) as writer:
                for i in table.entries:
                    writer.write(i)
        except exception.BedError:
            sys.stderr.write('Incorrect BED file.\n')


def rmout2bed_parser(subparsers):
    """
    Parser for the rmout2bed tool.
    """
    parser = subparsers.add_parser(
        'rmout2bed',
        help='convert a RepeatMasker out file to the BED format',
        description='Convert a RepeatMasker out file to the BED '
                    'format.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('rmout_file',
                        help='a RepeatMasker out file')
    parser.add_argument('bed_file',
                        help='the output BED file')

    # optional arguments
    parser.add_argument('-c', '--color', default='class',
                        choices=('class',
                                 'identity',
                                 'class_identity'),
                        help='how to choose colors of BED repeat '
                             'records')
    parser.add_argument('-n', '--name', default='id',
                        choices=('id',
                                 'name',
                                 'class',
                                 'family',
                                 'class_family'),
                        help='how to choose names of BED repeat '
                             'records')
    parser.add_argument('-s', '--short', action='store_true',
                        help='output only repeat loci (the output is '
                             'a BED3 file)')


def rmout2bed_launcher(args):
    """
    Launcher for the rmout2bed tool.
    """
    with open(args.rmout_file) as repeatmasker_file:
        with bed.Writer(args.bed_file) as bed_writer:
            rm_reader = repeatmasker.Reader(repeatmasker_file)
            for record in rm_reader.repeats():
                bed_record = repeatmasker.rmout2bed_record(
                    record, args.name, args.color, args.short)
                bed_writer.write(bed_record)


def gfftagstat_parser(subparsers):
    """
    Parser for the gfftagstat tool.
    """
    parser = subparsers.add_parser(
        'gfftagstat',
        help='get statistics on attribute tags of a GFF3 file',
        description='Get statistics on attribure tags of a GFF file.'
    )
    parser.add_argument('gff_file', help='a GFF3 file')

    # optional arguments
    parser.add_argument('-s', '--source', default=None,
                        help='filter GFF3 features by the specified '
                             'source')
    parser.add_argument('-t', '--type', default=None,
                        help='filter GFF3 features by the specified '
                             'type')


def gfftagstat_launcher(args):
    """
    Launcher for the gfftagstat tool.
    """
    with open(args.gff_file) as gff_file:
        tag_stats = gff3.analyze_tags(gff_file, args.source, args.type)
        is_filtered = args.source is not None or args.type is not None
        total_count = tag_stats['total']
        filtered_count = tag_stats['filtered']
        if is_filtered:
            print('Tag\tCount\t% Filtered\t% Total')
            output_template = '{0}\t{1}\t{2:.2f}\t{3:.2f}'
        else:
            print('Tag\tCount\t% Total')
            output_template = '{0}\t{1}\t{3:.2f}'

        print(output_template.format(
            '#records',
            total_count,
            100.0 * total_count / filtered_count,
            100
        ))
        if is_filtered:
            print(output_template.format(
                '#filtered',
                filtered_count,
                100,
                100.0 * filtered_count / total_count
            ))
        tag_counts = tag_stats['tag_counts']
        for tag in sorted(tag_counts.keys()):
            print(output_template.format(
                tag,
                tag_counts[tag],
                100.0 * tag_counts[tag] / filtered_count,
                100.0 * tag_counts[tag] / total_count
            ))


def gff2to3_parser(subparsers):
    """
    Parser for the gff2to3 tool.
    """
    parser = subparsers.add_parser(
        'gff2to3',
        help='convert a GFF2 file to the GFF3 format',
        description='Convert a GFF2 file to the GFF3 format'
    )
    parser.add_argument('gff2_file', help='a GFF2 file')
    parser.add_argument('output_file', help='the output GFF3 file')

    # optional arguments
    parser.add_argument('-i', '--ignore_incorrect_records',
                        action='store_true',
                        help='ignore incorrect records in the '
                             'specified input GFF2 file')


def gff2to3_launcher(args):
    """
    Launcher for the gff2to3 tool.
    """
    with open(args.gff2_file) as input_file:
        with open(args.output_file, 'w') as output_file:
            gff3.gff2to3(input_file, output_file,
                         not args.ignore_incorrect_records)


def snpeff2pph_parser(subparsers):
    """
    Parser for the snpeff2pph tool.
    """
    parser = subparsers.add_parser(
        'snpeff2pph',
        help='create a PolyPhen2 input file from an snpEff-annotated '
             'VCF file',
        description='Create an input file of aminoacid substitutions '
                    'for PolyPhen2 from a VCF file of variants '
                    'annotated with snpEff'
    )

    parser.add_argument('vcf_file', help='an snpEff-annotated VCF '
                                         'file')
    parser.add_argument('output_file', help='the output file in the '
                                            'PolyPhen2 format')


def snpeff2pph_launcher(args):
    """
    Launcher for the snpeff2pph tool.
    """
    with open(args.vcf_file) as vcf_file:
        with open(args.output_file, 'w') as output_file:
            template = '\t'.join(['{}'] * 4) + '\n'
            reader = vcf.Reader(vcf_file)
            for variant in reader:
                if 'ANN' in variant.INFO:
                    for effect in variant.INFO['ANN']:
                        snpeff_ann = snpeff.parse_snpeff_ann(effect)
                        if snpeff_ann.annotation == 'missense_variant':
                            output_file.write(template.format(
                                snpeff_ann.feature_id,
                                snpeff_ann.hgvs_p.pos,
                                snpeff.aa_code[snpeff_ann.hgvs_p.ref],
                                snpeff.aa_code[snpeff_ann.hgvs_p.alt]))


def gff2bed_parser(subparsers):
    """
    Parser for the gff2bed tool.
    """
    parser = subparsers.add_parser(
        'gff2bed',
        help='convert a GFF3 file to the BED format',
        description='Convert a GFF3 file to the BED format.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('gff_file', help='a GFF3 file')
    parser.add_argument('type', help='type of features to be '
                                     'processed')
    parser.add_argument('output_file',
                        help='the output file in the BED format')

    # optional arguments
    parser.add_argument('-a', '--attributes', nargs='*',
                        help='attributes to include to the output BED '
                             'file')
    parser.add_argument('-n', '--name_tag', default=None,
                        help='an attribute tag of a feature name')
    parser.add_argument('-m', '--missing_value', default='NA',
                        help='the missing tag value')
    parser.add_argument('-g', '--genes', action='store_true',
                        help='output a BED12 file of genes')
    parser.add_argument('-p', '--parent_tag', default='Parent',
                        help='an attribute tag of exon genes')
    parser.add_argument('--no_order_check', action='store_true',
                        help='do not check the order of GFF3 file '
                             'records')


def gff2bed_launcher(args):
    if args.genes:
        bed.convert_gff2bed_gene(args.gff_file, args.output_file,
                                 args.type, args.parent_tag,
                                 args.no_order_check)
    else:
        bed.convert_gff2bed(args.gff_file, args.output_file,
                            args.type, args.name_tag,
                            args.missing_value, args.attributes)


def snpeff2bed_parser(subparsers):
    """
    Parser for the snpeff2bed tool.
    """
    parser = subparsers.add_parser(
        'snpeff2bed',
        help='convert an snpEff-annotated VCF file to the BED format',
        description='Convert an snpEff-annotated VCF file to the BED '
                    'format with extra columns that describe variant '
                    'effects'
    )
    parser.add_argument('vcf_file', help='an snpEff-annotated VCF '
                                         'file')
    parser.add_argument('bed_file', help='the output BED file of '
                                         'annotated variants')
    parser.add_argument('--bed3', action='store_true',
                        help='convert to the BED3 format')


def snpeff2bed_launcher(args):
    """
    Launcher for the snpeff2bed tool.
    """
    snpeff.convert_snpeff2bed(args.vcf_file, args.bed_file, args.bed3)


def vcfgeno2bed_parser(subparsers):
    """
    Parser for the vcfgeno2bed tool.
    """
    parser = subparsers.add_parser(
        'vcfgeno2bed',
        help='extract genotype counts from a VCF file in the BED3+ '
             'format',
        description='Given a VCF file, extract genotype counts from '
                    'it and write them to the specified file in the '
                    'BED3+ format.'
    )
    parser.add_argument('vcf_file', help='a VCF file')
    parser.add_argument('output_file', help='the output BED3+ file of '
                        'genotype counts')
    parser.add_argument('-i', '--individuals', default=None,
                        help='a file with the list of individuals to '
                             'be considered for genotype counting')


def vcfgeno2bed_launcher(args):
    """
    Launcher for the vcfgeno2bed tool.
    """
    if args.individuals is not None:
        individuals = []
        with open(args.individuals) as individuals_file:
            for line in individuals_file:
                individuals.append(line.rstrip())
    else:
        individuals = None
    with open(args.vcf_file) as input_file:
        variants.convert_vcf2genotypes(input_file, args.output_file,
                                       individuals)


def vcfeffect2bed_parser(subparsers):
    """
    Parser for the vcfeffect2bed tool.
    """
    parser = subparsers.add_parser(
        'vcfeffect2bed',
        help='create a BED file of genotype effects from an '
             'snpEff-annotated VCF file',
        description='Given an snpEff-annotated VCF file, extract its '
                    'sample genotype effects.',
        epilog='Genotype values:\n'
               '\tREFHET - a heterozygote with one reference allele\n'
               '\tCOMHET - a heterozygote with both alternative '
               'alleles\n'
               '\tALTHOM - an alternative homozygote\n',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('vcf_file', help='an snpEff-annotated VCF '
                                         'file')
    parser.add_argument('output_file', help='the output BED3+ file of '
                        'sample effects')
    parser.add_argument('-i', '--impacts', nargs='+', choices=[
        'HIGH', 'MODERATE', 'LOW', 'MODIFIER'], default={'HIGH',
                                                         'MODERATE',
                                                         'LOW',
                                                         'MODIFIER'},
                        help='impacts of effects to be reported')
    parser.add_argument('-g', '--genotypes', nargs='+', choices=[
        'REFHET', 'COMHET', 'ALTHOM'], default={'REFHET', 'COMHET',
                                                'ALTHOM'})
    parser.add_argument('--ignore_errors', action='store_true',
                        help="ignore errors in an input file")


def vcfeffect2bed_launcher(args):
    """
    Launcher for the vcfeffect2bed tool.
    """
    try:
        snpeff.convert_vcfeffect2bed(args.vcf_file, args.output_file,
                                     impacts=args.impacts,
                                     genotypes=args.genotypes,
                                     ignore_errors=args.ignore_errors)
    except (KeyError, exception.BioformatsError):
        # let a user know if something goes wrong
        with open(args.output_file, "a") as output_file:
            output_file.write("ERROR!\n")


def flanknfilter_parser(subparsers):
    """
    Parser for the flanknfilter tool.
    """
    parser = subparsers.add_parser(
        "flanknfilter",
        help="filter features from a BED or VCF file by having N's "
             "in their flanking regions",
        description="Given features from a BED or VCF file, check if "
                    "they contain N's in their flanking regions of "
                    "the specified length.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", help="an input file of features "
                                           "to be filtered")
    parser.add_argument("fasta_file", help="a FASTA file of sequences "
                                           "the features are related "
                                           "to")
    parser.add_argument("output_file", help="an output file of "
                                            "filtered features")
    parser.add_argument("-t", "--type", choices=["bed", "vcf"],
                        default="bed", help="the input file type")
    parser.add_argument("-l", "--length", type=int, default=100,
                        help="the flanking region length")
    parser.add_argument("-s", "--strict", action="store_true",
                        help="require flanks to have exactly the "
                             "specified length (it may be shorter if "
                             "a feature is located near a sequence "
                             "start or end)")


def flanknfilter_launcher(args):
    """
    Launcher for the flanknfilter tool.
    """
    feature_filter = fasta.FlankNFilter(args.fasta_file,
                                        args.length)
    if args.type == 'bed':
        feature_filter.filter_bed(args.input_file, args.output_file,
                                  args.strict)
    else:
        feature_filter.filter_vcf(args.input_file, args.output_file,
                                  args.strict)


def interval2bed_parser(subparsers):
    """
    Parser for the interval2bed tool.
    """
    parser = subparsers.add_parser(
        'interval2bed',
        help='convert an interval file to the BED format',
        description='Convert a file in the interval format to the '
                    'BED format.'
    )
    parser.add_argument('interval_file', help='an interval file')
    parser.add_argument('bed_file', help='the output BED file')


def interval2bed_launcher(args):
    """
    Launcher for the interval2bed tool.
    """
    interval.convert_interval2bed(args.interval_file, args.bed_file)


def vcf2bed_parser(subparsers):
    """
    Parser for the vcf2bed tool.
    """
    parser = subparsers.add_parser(
        'vcf2bed',
        help='convert a VCF file to the BED format',
        description='Convert a file in the VCF format to the BED '
                    'format.'
    )
    parser.add_argument('vcf_file', help='a VCF file')
    parser.add_argument('bed_file', help='the output BED file')


def vcf2bed_launcher(args):
    """
    Launcher for the vcf2bed tool.
    """
    variants.vcf2bed(args.vcf_file, args.bed_file)
