#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
import gzip
import logging
from collections import namedtuple

logging.basicConfig()
logger = logging.getLogger(__name__)

record_names = ('chrom', 'pos', 'ref', 'alt', 'counts')

Record = namedtuple('Record', record_names)


class Reader(object):
    """
    This class provides routines to read multiple allele counts
    produced from VCFtools frequency count files (*.frq.count).
    """

    def __init__(self, filename, gzipped=False):
        """
        Create a Reader object from the specified file.

        :param filename: a name of a multiple allele count file
        :param gzipped: is the specified file gzipped or not
        :type filename: str
        :type gzipped: bool
        """
        self.__filename = filename
        self.__gzipped = gzipped
        self.__lineno = 0
        self.__line = ''

    def alleles(self):
        """
        Return the iterator to allele counts of the variants.

        :return: the allele count record iterator
        """
        file_handler = open if not self.__gzipped else gzip.open
        with file_handler(self.__filename) as count_file:
            reader = csv.reader(count_file, delimiter='\t')
            for line in reader:
                line[4] = tuple(map(int, line[4:]))
                yield Record(*line[:5])
