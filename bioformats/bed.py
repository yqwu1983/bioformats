#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
import logging
from collections import namedtuple
from . import autosql
from . import gff3
from .exception import BedError

logging.basicConfig()
logger = logging.getLogger(__name__)

bed_columns = ('seq', 'start', 'end', 'name', 'score', 'strand',
               'thick_start', 'thick_end', 'color', 'block_num',
               'block_sizes', 'block_starts', 'extra')

bed_numeric_fields = (1, 2, 4, 6, 7, 9)

Record = namedtuple('Record', bed_columns)

bed_autosql_fields = (
    autosql.TableEntry(type='string', num=None, name='chrom',
                       desc='Reference sequence chromosome or '
                            'scaffold'),
    autosql.TableEntry(type='uint', num=None, name='chromStart',
                       desc='Start position of feature on chromosome'),
    autosql.TableEntry(type='uint', num=None, name='chromEnd',
                       desc='End position of feature on chromosome'),
    autosql.TableEntry(type='string', num=None, name='name',
                       desc='Name of feature'),
    autosql.TableEntry(type='uint', num=None, name='score',
                       desc='Score'),
    autosql.TableEntry(type='char', num=1, name='strand',
                       desc='+ or - for strand'),
    autosql.TableEntry(type='uint', num=None, name='thickStart',
                       desc='Coding region start'),
    autosql.TableEntry(type='uint', num=None, name='thickEnd',
                       desc='Coding region end'),
    autosql.TableEntry(type='uint', num=None, name='reserved',
                       desc='Color set'),
    autosql.TableEntry(type='int', num=None, name='blockCount',
                       desc='The number of blocks in feature'),
    autosql.TableEntry(type='int[blockCount]', num=None,
                       name='blockSizes',
                       desc='Block sizes'),
    autosql.TableEntry(type='int[blockCount]', num=None,
                       name='chromStarts',
                       desc='Block start positions')
)


def is_score(x):
    """
    Given a value, check if it can represent a BED score.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent a BED score
    :rtype: bool
    """
    return autosql.is_int(x) and (0 <= int(x) <= 1000)


def is_strand(x):
    """
    Given a value, determine if it can represent a BED strand.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent a BED strand.
    :rtype: bool
    """
    return x in ('+', '-')


def is_coord(x):
    """
    Given a value, determine if it can represent a BED coordinate.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent a BED coordinate
    :rtype: bool
    """
    return autosql.is_int(x) and (int(x) >= 0)


def is_itemrgb(value):
    """
    Given a value, determine if it can represent a BED RGB color
    value.

    :param value: a value to check
    :type value: str
    :return: if the specified value can represent a BED RGB color value
    :rtype: bool
    """
    color_values = value.split(',', 2)
    if len(color_values) < 3:
        return False
    if not all([autosql.is_int(x) for x in color_values]):
        return False
    color_codes = [int(x) for x in color_values]
    for i in color_codes:
        if not (0 <= i <= 255):
            return False
    return True


def is_block_count(value):
    """
    Given a value, determine if it can represent a BED block count.

    :param value: a value to check
    :type value: str
    :return: if the specified value can represent a BED block count
    :rtype: bool
    """
    return autosql.is_int(value) and (int(value) > 0)


def is_block_sizes(value):
    """
    Given a value, determine if it can represent BED block sizes.

    :param value: a value to check
    :type value: str
    :return: if the specified value can represent BED block sizes
    :rtype: bool
    """
    size_values = value.split(',')
    if not all([autosql.is_int(x) for x in size_values]):
        return False
    size_numbers = [int(x) for x in size_values]
    for value in size_numbers:
        # block sizes must be positive numbers
        if not value > 0:
            return False
    return True


def is_block_starts(value):
    """
    Given a value, determine if it can represent BED block sizes.

    :param value: a value to check
    :type value: str
    :return: if the specified value can represent BED block sizes
    :rtype: bool
    """
    size_values = value.split(',')
    if not all([autosql.is_int(x) for x in size_values]):
        return False
    size_numbers = [int(x) for x in size_values]
    # the first block must start with 0
    if size_numbers[0] != 0:
        return False
    if len(size_numbers) > 1:
        # block starts must be ascending numbers
        for value, y in zip(size_numbers[:-1], size_numbers[1:]):
            if not (value < y):
                return False
    return True

bed_field_check = (
    lambda x: True,
    is_coord,
    is_coord,
    lambda x: True,
    is_score,
    is_strand,
    is_coord,
    is_coord,
    is_itemrgb,
    is_block_count,
    is_block_sizes,
    is_block_starts
)


class Reader(object):
    """
    This class implements a parser to read data from a file in the
    BED format.
    """

    def __init__(self, handle):
        """
        Given a handle of a file, create a BED reader object to read
        data from it.

        :param handle: a handle of a BED file
        """
        self.__reader = csv.reader(handle, delimiter='\t')
        self.__line_parts = []
        self.__bed_col = 12     # the number of BED columns
        self.__aux_col = 0      # the number of auxiliary columns

    def records(self, check_order=False):
        """
        Iterate through records in the BED file the object was
        created from.

        :param check_order: check if BED records are sorted; if they
            are not, then raise the exception
        :type check_order: bool
        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        prev_seq = ''
        prev_start = -1
        for self.__line_parts in self.__reader:
            new_record = self.__parse_bed_line()
            if check_order and ((prev_seq > new_record.seq) or (
                    prev_start > new_record.start)):
                logger.error('line %d: BED record order violated',
                             self.__reader.line_num)
                raise BedError
            prev_seq = new_record.seq
            prev_start = new_record.start
            yield new_record


    def __get_bed_format(self):
        """
        Determine the number of BED columns and the number of extra
        columns in a BED file.

        :return: a tuple of two numbers: the number of BED columns
            and the number of extra columns
        :rtype: tuple
        """
        i = 0
        for i, value in enumerate(self.__line_parts):
            if i < 12:
                if not bed_field_check[i](value):
                    i -= 1
                    break
            else:
                return 12, len(self.__line_parts) - 12
        return i + 1, len(self.__line_parts) - (i + 1)

    @property
    def bed_columns(self):
        """
        Get the number of BED columns in a BED file being read.

        :return: the number of BED columns
        :rtype: int
        """
        return self.__bed_col

    @property
    def aux_columns(self):
        """
        Get the number of auxiliary columns in a BED file being read.

        :return: the number of auxiliary BED columns
        :rtype: int
        """
        return self.__aux_col

    def __parse_bed_line(self):
        """
        Parse the current line from the BED file.

        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        bed_col, aux_col = self.__get_bed_format()
        self.__bed_col = min(self.__bed_col, bed_col)
        self.__aux_col = max(self.__aux_col, aux_col)
        if self.__bed_col == 7:
            # thickStart and thickEnd columns must be present together
            self.__bed_col = 6
            self.__aux_col += 1
        elif 10 <= self.bed_columns < 12:
            # blockCount, bloclSizes and blockStarts columns must be
            # present together
            self.__aux_col += self.__bed_col - 9
            self.__bed_col = 9

        if self.__bed_col < 3:
            # The first three columns of a BED file are mandatory.
            logger.debug('incorrect BED line %d',
                         self.__reader.line_num)
            raise BedError

        # convert numeric values: start, end, score, thick_start,
        # thick_end, block_num, blocks_sizes and block_starts; the
        # given BED line may contain the lesser number of columns,
        # so we adjust the tuple of numeric value positions
        line_numeric_pos = [x for x in bed_numeric_fields
                            if x < self.__bed_col]
        for i in line_numeric_pos:
            self.__line_parts[i] = int(self.__line_parts[i])

        # form the tuple to be returned as a result
        bed_parts = self.__line_parts[:self.__bed_col]
        bed_parts += ([None] * (12 - self.__bed_col))
        aux_parts = self.__line_parts[self.__bed_col:]
        result = Record(*bed_parts, extra=aux_parts)

        return result


class Writer(object):
    """
    The class implements writing to a file in the BED format.
    """
    def __init__(self, filename):
        """
        Given a name of a file, create a BED writer object to write
        data to it.

        :param filename: a name of a file to write BED records to
        :type filename: str
        """
        self.__filename = filename

    def __enter__(self):
        self.__output = open(self.__filename, 'w')
        return self

    def write(self, bed_record):
        """
        Given a BED record, write it to the file specified when the
        object was creted.

        :param bed_record: a BED record to be written to the file
        :type: Record
        """
        bed_record = [x for x in bed_record if x is not None]
        num_fields = len(bed_record)
        # check if the last column contains any values
        if not bed_record[-1]:
            num_fields -= 1
        else:
            num_fields += (len(bed_record[-1]) - 1)
            bed_record = bed_record[:-1] + bed_record[-1]
        template = '\t'.join(['{}'] * num_fields) + '\n'
        self.__output.write(template.format(*bed_record))

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()


def get_autosql_table(bed_reader, name='Table', desc='Description',
                      lines=None):
    """
    Given a BED reader, process entries from the associated BED file
    and choose appropriate autoSql types.

    :param bed_reader: a BED reader
    :param name: a table name
    :param desc: a table description
    :param lines: the number of lines to analyze
    :type bed_reader: Reader
    :type name: str
    :type desc: str
    :type lines: int
    :return: an autoSql table
    :rtype: autosql.Table
    """
    # process the first record from the BED reader to determine the
    # number of columns in the associated file
    first_record = next(bed_reader.records())
    num_aux_columns = bed_reader.aux_columns

    # create autoSql type classifiers for extra columns; for BED
    # columns, we already know their format and do not need classifiers
    column_types = []
    for i in range(num_aux_columns):
        new_classifier = autosql.Classifier()
        new_classifier.add_value(str(first_record.extra[i]))
        column_types.append(new_classifier)

    # we have already checked the first line
    if lines is not None:
        lines -= 1
    # process entries from the specified BED reader
    for line_num, record in enumerate(bed_reader.records()):
        if lines is not None and not (line_num < lines):
            break
        for i in range(num_aux_columns):
            column_types[i].add_value(str(record.extra[i]))

    # form the table of autoSql entries
    entries = []
    for i in range(bed_reader.bed_columns):
        entries.append(bed_autosql_fields[i])

    for i in range(bed_reader.aux_columns):
        entries.append(autosql.TableEntry(
            type=column_types[i].data_type,
            name='column_{}'.format(i+1),
            desc='Column #{} with {} values'.format(
                i+1, column_types[i].data_type),
            num=None
        ))

    table_scheme = autosql.Table(name=name, desc=desc, entries=entries)
    return table_scheme


def get_blocks(start, end):
    """
    Given start and end coordinates of regions, return their starts
    and sizes in the BED block format.

    :param start: start positions of regions
    :param end: end positions of regions
    :return: a dictionary of two lists: block starts and sizes
    :rtype: dict
    """
    result = ([], [])

    if len(start) != len(end):
        logger.error('unequal number of start and end coordinates %d '
                     'and %d', len(start), len(end))
        raise BedError
    feature_start = start[0]
    for i, j in zip(start, end):
        result[0].append(i - feature_start)
        result[1].append(j - i + 1)

    return result


def convert_gff2bed_gene(gff3_file, bed_file, exon_type='exon',
                         parent_tag='Parent', ignore_order=False):
    """
    Convert a specified GFF3 file of gene exons to the BED12 format
    considering gene structure.

    :param gff3_file: a name of an input GFF3 file of gene exons
    :param bed_file: a name of the output BED12 file
    :param exon_type: a type of GFF3 exon records
    :param parent_tag: a tag of GFF3 exon records encoding the gene
        they belong to
    :param ignore_order: do not check the order of GFF3 records
    :type bed_file: str
    :type gff3_file: str
    :type exon_type: str
    :type parent_tag: str
    :type ignore_order: bool
    """
    total_exons = 0
    total_genes = 0
    with open(gff3_file) as input_file:
        gff_reader = gff3.Reader(input_file)
        with Writer(bed_file) as bed_writer:
            # process the first exon
            gff_record_iterator = gff_reader.records(
                    check_order=not ignore_order)
            exon = next(gff_record_iterator)
            while exon.type != exon_type:
                exon = next(gff_record_iterator)
            cur_gene_seq = exon.seqid
            cur_gene = exon.attributes[parent_tag]
            cur_gene_start = exon.start
            cur_gene_end = exon.end
            exon_starts = [exon.start]
            exon_ends = [exon.end]
            total_exons += 1
            # start iterating through exon records of a GFF3 file
            for exon in gff_record_iterator:
                if exon.type != exon_type:
                    continue
                total_exons += 1
                if exon.attributes[parent_tag] != cur_gene:
                    total_genes += 1
                    # we have read a gene, write it to output
                    block_starts, block_sizes = get_blocks(
                        exon_starts, exon_ends)
                    bed_record = Record(
                        seq=cur_gene_seq,
                        start=cur_gene_start - 1,
                        end=cur_gene_end,
                        name=cur_gene,
                        score=1000,
                        strand=exon.strand,
                        thick_start=cur_gene_start,
                        thick_end=cur_gene_end,
                        color='255,255,255',
                        block_num=len(exon_starts),
                        block_sizes=','.join(map(str, block_sizes)),
                        block_starts=','.join(map(str, block_starts)),
                        extra=[]
                    )
                    bed_writer.write(bed_record)
                    # update the exon data
                    cur_gene_seq = exon.seqid
                    cur_gene = exon.attributes[parent_tag]
                    cur_gene_start = exon.start
                    cur_gene_end = exon.end
                    exon_starts = [exon.start]
                    exon_ends = [exon.end]
                else:
                    # we have read another exon of the gene being
                    # processed
                    exon_starts.append(exon.start)
                    exon_ends.append(exon.end)
                    cur_gene_end = exon.end
            # process the last gene
            total_genes += 1
            block_starts, block_sizes = get_blocks(
                exon_starts, exon_ends)
            bed_record = Record(
                seq=cur_gene_seq,
                start=cur_gene_start - 1,
                end=cur_gene_end,
                name=cur_gene,
                score=1000,
                strand=exon.strand,
                thick_start=cur_gene_start,
                thick_end=cur_gene_end,
                color='255,255,255',
                block_num=len(exon_starts),
                block_sizes=','.join(map(str, block_sizes)),
                block_starts=','.join(map(str, block_starts)),
                extra=[]
            )
            bed_writer.write(bed_record)
    logger.info('%d exon records of %d genes processed', total_exons,
                total_genes)


def convert_gff2bed(gff3_file, bed_file, feature_type, name_tag=None,
                    missing_value='NA', attributes=None):
    """
    Convert a specified GFF3 file to the BED format considering its
    attributes.

    :param gff3_file: a name of an input GFF3 file
    :param bed_file: a name of the output BED file
    :param feature_type: a type of features to be processed
    :param name_tag: a tag which value to use as a feature name
    :param missing_value: the value to denote a missing attribute
    :param attributes: a list of attributes to be added to output
    :type gff3_file: str
    :type bed_file: str
    :type feature_type: str
    :type name_tag: str
    :type missing_value: str
    :type attributes: list
    """
    total_processed = 0
    total_bed = 0
    if attributes is None:
        attributes = []
    with open(gff3_file) as input_file:
        gff_reader = gff3.Reader(input_file)
        with Writer(bed_file) as bed_writer:
            for feature in gff_reader.records():
                total_processed += 1
                if feature.type != feature_type:
                    continue
                if name_tag is not None:
                    if name_tag in feature.attributes:
                        feature_name = feature.attributes[name_tag]
                    else:
                        feature_name = missing_value
                else:
                    feature_name = 'feature_{}'.format(total_bed)
                # prepare extra BED columns
                feature_extra = []
                for i in attributes:
                    if i in feature.attributes:
                        feature_extra.append(feature.attributes[i])
                    else:
                        feature_extra.append('NA')
                # prepare a BED record
                bed_record = Record(
                    seq=feature.seqid,
                    start=feature.start - 1,
                    end=feature.end,
                    name=feature_name,
                    score=1000,
                    strand=feature.strand,
                    thick_start=feature.start - 1,
                    thick_end=feature.end,
                    color=None,
                    block_num=None,
                    block_sizes=None,
                    block_starts=None,
                    extra=feature_extra
                )
                bed_writer.write(bed_record)
                total_bed += 1
    logger.info('%d BED records of %d GFF3 records processed',
                total_bed, total_processed)
