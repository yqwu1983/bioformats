#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from builtins import range  # pylint:disable=redefined-builtin
from collections import namedtuple
from .exception import BlastTabError
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)


class BlastTab(object):
    """
    This class implements a parser to read alignment files in the
    BLAST tabular format.
    """

    blast_field_names = ('query', 'subject', 'identity', 'length',
                         'mismatches', 'gap_openings',
                         'q_start', 'q_end', 's_start', 's_end',
                         'e_value', 'bit_score')

    Alignment = namedtuple('Alignment', blast_field_names)

    def __init__(self, filename):
        """
        Given a name of a file, create a BLAST tabular format parser
        object to read data from it.

        :param filename: a name of a file in the BLAST tabular format
        :type filename: str
        """
        self.__filename = filename
        self.__lineno = 0
        self.__line = None

    def alignments(self):
        """
        Iterate through alignments in the file the object was created
        from.

        :return: an alignment for the file the object was created from
        :rtype: Blast.Alignment
        """
        with open(self.__filename) as blast_file:
            for self.__line in blast_file:
                self.__lineno += 1
                if not self.__line.startswith('#'):
                    # if the line starts with '#', then it is a
                    # comment and we skip it
                    yield self.__parse_blast_line()

    def __parse_blast_line(self):
        """
        Parse the current line from the BLAST tabular file.

        :return: an alignment from the file the object was created from
        :rtype: Blast.Alignment
        """
        line_parts = self.__line.split('\t', 12)

        # check if the line contains the proper number of columns
        if len(line_parts) < 12:
            logging.error('line %d: the incorrect number of '
                          'columns', self.__lineno)
            raise BlastTabError

        # convert numeric values of identity, e-value and bit score
        # to float numbers
        for i in (2, 10, 11):
            try:
                line_parts[i] = float(line_parts[i])
            except ValueError:
                logging.error(
                    'line %d: the incorrect numerical value '
                    '%s', self.__lineno, line_parts[i])
                raise BlastTabError

        # convert numeric values of alignment length, the number of
        # mismatches, the number of gap openings and query and subject
        # coordinates
        for i in range(3, 10):
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logging.error(
                    'line %d: the incorrect integer value '
                    '%s', self.__lineno, line_parts[i])
                raise BlastTabError

        return BlastTab.Alignment(*line_parts)