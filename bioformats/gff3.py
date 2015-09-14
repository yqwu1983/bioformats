#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
import logging
from collections import defaultdict, namedtuple, OrderedDict
from future.utils import iteritems
from .exception import Gff3Error

logging.basicConfig()
logger = logging.getLogger(__name__)

gff3_columns = ('seqid', 'source', 'type', 'start', 'end', 'score',
                'strand', 'phase', 'attributes', )

Record = namedtuple('Record', gff3_columns)


class Reader(object):
    """
    This class implements a parser to read data from a file in the
    GFF3 format.
    """
    def __init__(self, handle):
        """
        Given a handle of a file, create a GFF3 reader object to read
        data from it.

        :param handle: a handle of a GFF3 file
        """
        self.__reader = csv.reader(handle, delimiter='\t')
        self.__line_parts = []

    def records(self):
        """
        Iterate through records in the GFF3 file the object was
        created from.

        :return: a record from the GFF3 file the object was created
        from
        :rtype: Record
        """
        # skip the first line of the file
        next(self.__reader)
        for self.__line_parts in self.__reader:
            yield self.__parse_gff3_line()

    def __parse_gff3_line(self):
        """
        Parse the current line from the GFF3 file.

        :return: a record from the GFF3 file the object was created
            from
        :rtype: Record
        """
        if len(self.__line_parts) < 8:
            logger.error('line %d: the incorrect number of columns - '
                         '%d', self.__reader.line_num,
                         len(self.__line_parts))
            raise Gff3Error

        if len(self.__line_parts) > 9:
            # in the attributes column, some values may contain tab
            # characters that leads to multiple fields; so we
            # concatenate these fields to a single field
            self.__line_parts[8] = '\t'.join(self.__line_parts[8:])
            self.__line_parts = self.__line_parts[:9]

        # convert numeric values: start, end, score (if present) and
        # phase (if present)
        line_int_pos = [3, 4]
        if self.__line_parts[7] != '.':
            line_int_pos += [7]
        for i in line_int_pos:
            try:
                self.__line_parts[i] = int(self.__line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect numeric value '
                             '%s', self.__reader.line_num,
                             self.__line_parts[i])
                raise Gff3Error
        if self.__line_parts[5] != '.':
            # if the score is specified, try to convert it to a float
            # number value
            try:
                self.__line_parts[5] = float(self.__line_parts[5])
            except ValueError:
                logger.error('line %d: the incorrect float number '
                             'value %s', self.__reader.line_num,
                             self.__line_parts[5])
                raise Gff3Error

        if len(self.__line_parts) == 9:
            # parse the attributes
            attributes = OrderedDict()
            for x in self.__line_parts[8].split(';'):
                try:
                    tag, value = x.split('=', 2)
                except ValueError:
                    logger.error('line %d: the incorrect GFF3 '
                                 'attribute %s',
                                 self.__reader.line_num, x)
                    raise Gff3Error
                if not tag or not value:
                    logger.error('line %d: the incorrect GFF3 '
                                 'attribute %s',
                                 self.__reader.line_num, x)
                    raise Gff3Error
                attributes[tag] = value
            self.__line_parts[8] = attributes
        else:
            self.__line_parts += [None]

        return Record(*self.__line_parts)


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
        :type gff3_record: Record
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


def analyze_tags(handle, feature_source=None, feature_type=None):
    """
    Given a handle of a GFF3 file, collect statistics on attribute
    tags in it.

    :param handle: a handle of a GFF3 file
    :param feature_source: if specified, then only features of the
        specified source are considered
    :param feature_type: if specified, then only features of the
        specified type are considered
    :type feature_source: str
    :type feature_type: str
    :return: a dictionary of three values: the total number of
        processed features, the number of filtered features and the
        dictionary of attribute tag counts
    :rtype: dict
    """
    total_features = filtered_features = attr_features =  0
    tag_counts = defaultdict(int)
    for record in Reader(handle).records():
        total_features += 1
        if feature_source is not None and record.source != \
                feature_source:
            continue
        if feature_type is not None and record.type != feature_type:
            continue
        filtered_features += 1
        if record.attributes is not None:
            for tag in record.attributes:
                tag_counts[tag] += 1

    return {'total': total_features,
            'filtered': filtered_features,
            'tag_counts': tag_counts}
