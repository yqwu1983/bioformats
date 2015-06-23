#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
from collections import namedtuple
from future.utils import iteritems
from .exception import Gff3Error

logging.basicConfig()
logger = logging.getLogger(__name__)

gff3_columns = ('seqid', 'source', 'type', 'start', 'end', 'score',
                'strand', 'phase', 'attributes', )

Gff3Record = namedtuple('GffRecord', gff3_columns)


class Reader(object):
    """
    This class implements a parser to read data from a file in the
    GFF3 format.
    """
    def __init__(self, filename):
        """
        Given a name of a file, create a GFF3 reader object to read
        data from it.

        :param filename: a name of a GFF3 file
        :type filename: str
        """
        self.__filename = filename
        self.__lineno = 0
        self.__line = ''

    def records(self):
        """
        Iterate through records in the GFF3 file the object was
        created from.

        :return: a record from the GFF3 file the object was created
        from
        :rtype: Gff3Record
        """
        with open(self.__filename) as gff3_file:
            # check the first line of the file
            self.__line = next(gff3_file)
            self.__line = self.__line.rstrip()
            self.__lineno += 1
            if self.__line != '##gff-version 3':
                logger.error('line %d: the incorrect GFF3 header '
                             '%s'.format(self.__lineno, self.__line))
                raise Gff3Error
            for self.__line in gff3_file:
                self.__lineno += 1
                yield self.__parse_gff3_line()

    def __parse_gff3_line(self):
        """
        Parse the current line from the GFF3 file.

        :return: a record from the GFF3 file the object was created
            from
        :rtype: Gff3Record
        """
        line_parts = self.__line.rstrip().split('\t', 8)
        if len(line_parts) < 8:
            logger.error('line %d: the incorrect number of columns - '
                         '%d', self.__lineno, len(line_parts))
            raise Gff3Error

        # convert numeric values: start, end, score (if present) and
        # phase (if present)
        line_int_pos = [3, 4]
        if line_parts[7] != '.':
            line_int_pos += [7]
        for i in line_int_pos:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect numeric value '
                             '%s', self.__lineno, line_parts[i])
                raise Gff3Error
        if line_parts[5] != '.':
            # if the score is specified, try to convert it to a float
            # number value
            try:
                line_parts[5] = float(line_parts[5])
            except ValueError:
                logger.error('line %d: the incorrect float number '
                             'value %s', self.__lineno, line_parts[5])
                raise Gff3Error

        if len(line_parts) == 9:
            # parse the attributes
            attributes = dict()
            for x in line_parts[8].split(';'):
                tag, value = x.split('=', 2)
                attributes[tag] = value
            line_parts[8] = attributes
        else:
            line_parts += [None]

        return Gff3Record(*line_parts)


class Writer(object):
    """
    The class implements writing to a file in the GFF3 format.
    """
    def __init__(self, filename):
        """
        Given a name of a file, create a GFF3 writer object to write
        data to it.

        :param filename: a name of a file to write GFF3 records to
        :type filename: str
        """
        self.__filename = filename

    def __enter__(self):
        self.__output = open(self.__filename, 'w')
        self.__output.write('##gff-version 3\n')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()

    def write(self, gff3_record):
        """
        Given a GFF3 record, write it to the file specified when the
        object was created.

        :param gff3_record: a GFF3 record to be written to the file
        :return: Gff3Record
        """
        gff3_record = list(gff3_record)
        # prepare the attributes line
        if gff3_record[-1]:
            attributes = ['{}={}'.format(x, y) for x, y in
                          iteritems(gff3_record[-1])]
            gff3_record[-1] = ';'.join(attributes)
        else:
            gff3_record = gff3_record[:-1]

        template = '\t'.join(['{}'] * len(gff3_record)) + '\n'
        self.__output.write(template.format(*gff3_record))
