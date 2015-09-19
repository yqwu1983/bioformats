#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
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

    def __init__(self, handle):
        """
        Create a Reader object from the specified file.

        :param handle: a handle of a multiple allele count file
        """
        self.__handle = handle
        self.__lineno = 0
        self.__line = ''

    def alleles(self):
        """
        Return the iterator to allele counts of the variants.

        :return: the allele count record iterator
        """
        reader = csv.reader(self.__handle, delimiter='\t')
        for line in reader:
            line[4] = tuple(map(int, line[4:]))
            yield Record(*line[:5])
