#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple
from bioformats.exception import FrqCountReaderError
import gzip
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)


class Reader(object):
    """
    This class provides routines to process allele frequency files
    (*.frq.count) produced by VCFtools using its --counts option.
    """

    numeric_fields = (1, 2, 3)

    frq_count_names = ('chrom', 'pos', 'n_alleles', 'n_chr', 'ref',
                       'alt')

    Record = namedtuple('FrqCountLine', frq_count_names)

    def __init__(self, filename, gzipped=False):
        """
        Create an Reader object from the specified file.

        :param filename: a name of a VCFtools frequency count file
        :param gzipped: is the specified file gzipped or not
        :type filename: str
        :type gzipped: bool
        """
        self.__filename = filename
        self.__gzipped = gzipped
        self.__lineno = 0
        self.__line = ''

    def __parse_frq_count_line(self):
        """
        Parse a current line from a frequency count file.

        :return: a named tuple representing a frequecy count record
        :rtype: Reader.Record
        """
        line_parts = self.__line.split(None, 5)

        # check if the line contains all required fields
        if len(line_parts) < 6:
            logger.error('line %d: the incorrect number of '
                         'fields', self.__lineno)
            raise FrqCountReaderError

        # convert numeric values
        for i in Reader.numeric_fields:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect numeric value '
                             '%s', self.__lineno, line_parts[i])
                raise FrqCountReaderError

        # parse the parts containing allele counts
        for i in (4, 5):
            line_parts[i] = self.__parse_allele_counts(
                line_parts[i])

        return Reader.Record(*line_parts)

    def __parse_allele_counts(self, counts):
        """
        Parse allele counts from a current frequency count record

        :param counts: a line of allele counts
        :type counts: str
        :return: a dictionary which keys are alleles and values are
            their counts
        :rtype: dict
        """
        result = dict()

        for allele_record in counts.split():
            try:
                allele, count = allele_record.rsplit(':', 1)
            except ValueError:
                logger.error('line %d: the incorrect allele record '
                             '%s', self.__lineno, allele_record)
                raise FrqCountReaderError
            # check if the count value is a valid integer
            try:
                count = int(count)
            except ValueError:
                logger.error('line %d: the incorrect allele count '
                             'value %s', self.__lineno, count)
                raise FrqCountReaderError
            # check if the count allele value is unique, otherwise
            # raise the exception
            if allele not in result:
                result[allele] = count
            else:
                logger.error('line %d: multiple allele counts for '
                             'the allele %s', self.__lineno, allele)
                raise FrqCountReaderError

        return result

    def variants(self):
        """
        Return the iterator to allele counts of the variants.

        :return: the allele count record iterator
        """
        file_handler = open if not self.__gzipped else gzip.open
        with file_handler(self.__filename) as count_file:
            self.__lineno = 0
            for self.__line in count_file:
                self.__line = self.__line.rstrip()
                self.__lineno += 1
                if self.__line.startswith('CHROM'):
                    # skip the comment line
                    continue
                yield self.__parse_frq_count_line()