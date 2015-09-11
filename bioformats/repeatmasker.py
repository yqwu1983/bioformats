#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
import logging
from collections import namedtuple
from .exception import RepeatMaskerReaderError

logging.basicConfig()
logger = logging.getLogger(__name__)

repeatmasker_out_columns = ('sw_score', 'subst_perc', 'del_perc',
                            'ins_perc', 'query', 'query_start',
                            'query_end', 'query_past',
                            'is_complement', 'repeat_name',
                            'repeat_class', 'repeat_prior',
                            'repeat_start', 'repeat_end',
                            'is_overlapping')

Record = namedtuple('Record', repeatmasker_out_columns)


class Reader(object):
    """
    This class implements routines to read data from a RepeatMasker
    out file.
    """

    int_fields = (0, 5, 6)
    float_fields = (1, 2, 3)

    def __init__(self, handle):
        """
        Create a RepeatMasker out reader from the specified handle.

        :param handle: a handle of a RepeatMasker out file
        """
        self.__reader = csv.reader(handle, delimiter='\t')
        self.__line_parts = []

    def repeats(self):
        """
        Return an iterator to itetate through repeat records from the
        RepeatMasker out file.

        :return: an iterator to iterate through repeat records
        """
        # skip the first three lines of the input file because they
        # contain a header
        for _ in range(3):
            next(self.__reader)
        for self.__line_parts in self.__reader:
            yield self.__parse_rmout_line()

    def __parse_rmout_line(self):
        """
        Parse the curent line from the RepeatMasker out file.

        :return: a record from the RepeatMasker out file the object
            was created from
        :rtype: Record
        """
        if not (14 <= len(self.__line_parts) <= 15):
            logger.error('line %d: the incorrect number of '
                         'columns', self.__reader.line_num)
            raise RepeatMaskerReaderError

        # add the last column value 'None' if it is missing
        if len(self.__line_parts) == 14:
            self.__line_parts += [None]

        for i in self.int_fields:
            try:
                self.__line_parts[i] = int(self.__line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect integer value '
                             '%s', self.__reader.line_num,
                             self.__line_parts[i])
                raise RepeatMaskerReaderError

        for i in self.float_fields:
            try:
                self.__line_parts[i] = float(self.__line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect float value %s',
                             self.__reader.line_num,
                             self.__line_parts[i])
                raise RepeatMaskerReaderError

        return Record(*self.__line_parts)
