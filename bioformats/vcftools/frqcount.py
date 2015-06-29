#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from future.utils import iteritems
from collections import namedtuple, OrderedDict
from bioformats.exception import FrqCountReaderError
import gzip
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)

numeric_fields = (1, 2, 3)

frqcount_names = ('chrom', 'pos', 'n_alleles', 'n_chr', 'ref', 'alt')

FrqCountRecord = namedtuple('FrqCountRecord', frqcount_names)


class Reader(object):
    """
    This class provides routines to process allele frequency files
    (*.frq.count) produced by VCFtools using its --counts option.
    """

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
        for i in numeric_fields:
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

        return FrqCountRecord(*line_parts)

    def __parse_allele_counts(self, counts):
        """
        Parse allele counts from a current frequency count record

        :param counts: a line of allele counts
        :type counts: str
        :return: a dictionary which keys are alleles and values are
            their counts
        :rtype: dict
        """
        result = OrderedDict()

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
        with file_handler(self.__filename, 'rt') as count_file:
            self.__lineno = 0
            for self.__line in count_file:
                self.__line = self.__line.rstrip()
                self.__lineno += 1
                if self.__line.startswith('CHROM'):
                    # skip the comment line
                    continue
                yield self.__parse_frq_count_line()


class Writer(object):
    """
    The class implements writing to a file in the VCFtools allele
    frequency format.
    """
    def __init__(self, filename, gzipped=False):
        """
        Given a name of a file, create a VCFtools frequencty count
        writer object to write data to it.

        :param filename: a name of a file to write allele frequency
            counts to
        :param gzipped: produce the gzipped output file or not
        :type filename: str
        :type gzipped: bool
        """
        self.__filename = filename
        self.__gzipped = gzipped

    def __enter__(self):
        if self.__gzipped:
            self.__output = gzip.open(self.__filename, 'wt')
        else:
            self.__output = open(self.__filename, 'w')
        self.__output.write(
            'CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()

    def write(self, frqcount_record):
        """
        Given an allele frequency count record, write it to the file
        specified when the object was created.

        :param frqcount_record: an allele frequency count record to
            be written to the file
        :type frqcount_record: FrqCountRecord
        """
        # prepare alelle frequencies
        ref_freq = ['{}:{}'.format(x, y)
                    for x, y in iteritems(frqcount_record.ref)]
        alt_freq = ['{}:{}'.format(x, y)
                    for x, y in iteritems(frqcount_record.alt)]
        ref_freq = ref_freq[0]
        alt_freq = '\t'.join(alt_freq)
        # prepare the line to be written
        frqcount_record = list(frqcount_record)
        frqcount_record[4] = ref_freq
        frqcount_record[5] = alt_freq

        template = '\t'.join(['{}'] * len(frqcount_record)) + '\n'
        self.__output.write(template.format(*frqcount_record))
