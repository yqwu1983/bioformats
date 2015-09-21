#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
import logging
from collections import namedtuple
from . import autosql
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
                       name='blockStarts',
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

    def records(self):
        """
        Iterate through records in the BED file the object was
        created from.

        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        for self.__line_parts in self.__reader:
            yield self.__parse_bed_line()

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


def get_autosql_table(bed_reader, name='Table', desc='Description'):
    """
    Given a BED reader, process entries from the associated BED file
    and choose appropriate autoSql types.

    :param bed_reader: a BED reader
    :param name: a table name
    :param desc: a table description
    :type bed_reader: Reader
    :type name: str
    :type desc: str
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

    # process entries from the specified BED reader
    for record in bed_reader.records():
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
