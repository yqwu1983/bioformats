#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
from collections import namedtuple
from .exception import RepeatMaskerError

logging.basicConfig()
logger = logging.getLogger(__name__)

repeatmasker_out_columns = ('sw_score', 'subst_perc', 'del_perc',
                            'ins_perc', 'query', 'query_start',
                            'query_end', 'query_past',
                            'is_complement', 'repeat_name',
                            'repeat_class', 'repeat_prior',
                            'repeat_start', 'repeat_end', 'id',
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
        self.__handle = handle
        self.__line = ''
        self.__lineno = 0

    def repeats(self):
        """
        Return an iterator to itetate through repeat records from the
        RepeatMasker out file.

        :return: an iterator to iterate through repeat records
        """
        # skip the first three lines of the input file because they
        # contain a header
        for _ in range(3):
            next(self.__handle)

        for self.__line in self.__handle:
            self.__lineno += 1
            yield self.__parse_rmout_line()

    def __parse_rmout_line(self):
        """
        Parse the curent line from the RepeatMasker out file.

        :return: a record from the RepeatMasker out file the object
            was created from
        :rtype: Record
        """
        line_parts = self.__line.rstrip().split(None, 14)
        if len(line_parts) < 15:
            logger.error('line %d: the incorrect number of '
                         'columns', self.__lineno)
            raise RepeatMaskerError

        # add the last column value 'None' if it is missing
        if len(line_parts) == 15:
            line_parts += [None]

        for i in self.int_fields:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect integer value '
                             '%s', self.__lineno, line_parts)
                raise RepeatMaskerError

        for i in self.float_fields:
            try:
                line_parts[i] = float(line_parts[i])
            except ValueError:
                logger.error('line %d: the incorrect float value %s',
                             self.__lineno, line_parts[i])
                raise RepeatMaskerError

        return Record(*line_parts)
