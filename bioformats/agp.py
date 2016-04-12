#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import csv
import logging
from .exception import AgpError

logging.basicConfig()
logger = logging.getLogger(__name__)

comp_numeric_fields = (1, 2, 3, 6, 7)
gap_numeric_fields = (1, 2, 3, 5)

component_types = ('A', 'D', 'F', 'G', 'O', 'P', 'W', 'N', 'U')
gap_types = ('N', 'U')


class Reader(object):
    """
    Parser to read data from an AGP file.
    """

    def __init__(self, handle):
        """
        Given a handle of a file, create an AGP reader object to read
        records from it.

        :param handle: a handle of an AGP file
        """
        self.__reader = csv.reader(handle, delimiter="\t", )
        self.__line_parts = []

    def records(self):
        """
        Iterate through records in the AGP file the object was
        created from.

        :return: a record from the AGP file
        :rtype: tuple
        """
        for self.__line_parts in self.__reader:
            if not self.__line_parts[0].startswith('#'):
                agp_record = self.__parse_agp_record()
                yield agp_record

    def __parse_agp_record(self):
        """
        Parse a line from an AGP file.

        :return: a record from an AGP file
        :rtype: tuple
        """
        line_parts = self.__line_parts
        line_num = self.__reader.line_num
        if len(line_parts) != 9:
            logger.error('line %d: incorrect number of columns: %d '
                         'instead of 9', line_num, len(line_parts))
            raise AgpError

        if line_parts[4] not in component_types:
            logger.error('line %d: incorrect component type %s',
                         line_num, line_parts[4])
            raise AgpError

        if line_parts[4] in gap_types:
            numeric_fields = gap_numeric_fields
        else:
            numeric_fields = comp_numeric_fields

        for i in numeric_fields:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line %d: incorrect numeric value '
                             '%s', line_num, line_parts[i])
                raise AgpError

        return line_parts
