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
        self.__handle = handle
        self.__lineno = 0
        self.__line_parts = []

    def records(self):
        """
        Iterate through records in the BED file the object was
        created from.

        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        reader = csv.reader(self.__handle, delimiter='\t')
        for self.__line_parts in reader:
            self.__lineno += 1
            yield self.__parse_bed_line()

    def __parse_bed_line(self):
        """
        Parse the current line from the BED file.

        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        line_parts = self.__line_parts

        # convert numeric values: start, end, score, thick_start,
        # thick_end, block_num, blocks_sizes and block_starts; the
        # given BED line may contain the lesser number of columns,
        # so we adjust the tuple of numeric value positions
        line_numeric_pos = [x for x in bed_numeric_fields
                            if x < len(line_parts)]
        for i in line_numeric_pos:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect numeric value '
                             '%s', self.__lineno, line_parts[i])
                raise BedError

        # form the tuple to be returned as a result
        line_parts += ([None] * (len(bed_columns) - len(line_parts)))
        result = Record(*line_parts)

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
        num_fields = len(bed_record)
        # check if the last column contains any values
        if bed_record.extra is None:
            num_fields -= 1
        template = '\t'.join(['{}'] * num_fields) + '\n'
        self.__output.write(template.format(*bed_record))

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()


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


def is_itemrgb(x):
    """
    Given a value, determine if it can represent a BED RGB color
    value.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent a BED RGB color value
    :rtype: bool
    """
    color_values = x.split(',', 2)
    if len(color_values) < 3:
        return False
    if not all([autosql.is_int(x) for x in color_values]):
        return False
    color_codes = [int(x) for x in color_values]
    for i in color_codes:
        if not (0 <= i <= 255):
            return False
    return True


def is_block_count(x):
    """
    Given a value, determine if it can represent a BED block count.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent a BED block count
    :rtype: bool
    """
    return autosql.is_int(x) and (int(x) > 0)


def is_block_sizes(x):
    """
    Given a value, determine if it can represent BED block sizes.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent BED block sizes
    :rtype: bool
    """
    size_values = x.split(',')
    if not all([autosql.is_int(x) for x in size_values]):
        return False
    size_numbers = [int(x) for x in size_values]
    for x in size_numbers:
        # block sizes must be positive numbers
        if not x > 0:
            return False
    return True


def is_block_starts(x):
    """
    Given a value, determine if it can represent BED block sizes.

    :param x: a value to check
    :type x: str
    :return: if the specified value can represent BED block sizes
    :rtype: bool
    """
    size_values = x.split(',')
    if not all([autosql.is_int(x) for x in size_values]):
        return False
    size_numbers = [int(x) for x in size_values]
    # the first block must start with 0
    if size_numbers[0] != 0:
        return False
    if len(size_numbers) > 1:
        # block starts must be ascending numbers
        for x, y in zip(size_numbers[:-1], size_numbers[1:]):
            if not (x < y):
                return False
    return True
