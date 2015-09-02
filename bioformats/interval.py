#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple
from .exception import IntervalError
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)

interval_fields = ('seq', 'start', 'end')
Record = namedtuple('Record', interval_fields)


class Reader(object):
    """
    This class implements a parser to read data from a file in the
    interval data format.
    """
    def __init__(self, filename):
        """
        Given a name of a file, create an interval format reader
        object to read data from it.

        :param filename: a name of an interval format file
        :type filename: str
        """
        self.__filename = filename
        self.__lineno = 0
        self.__line = ''
        self.__seq = ''

    def intervals(self):
        """
        Iterate through records in the interval format file the
        object was created from

        :return: a record from the interval format file the object
            was created from
        :rtype: Record
        """
        with open(self.__filename) as interval_file:
            for self.__line in interval_file:
                self.__line = self.__line.rstrip()
                self.__lineno += 1
                if self.__line.startswith('>'):
                    self.__seq = self.__line[1:]
                elif self.__line:
                    # parse the interval line
                    try:
                        start, end = (int(x) for x in
                                      self.__line.split('-', 2))
                    except (TypeError, ValueError):
                        logger.error('line %d: an incorrect interval '
                                     'line %s', self.__lineno,
                                     self.__line)
                        raise IntervalError
                    yield Record(self.__seq, start, end)


class Writer(object):
    """
    The class implements writing to a file in the interval format.
    """
    def __init__(self, filename):
        """
        Given a name of a file, create an interval format writer
        object to write data from it.

        :param filename: a name of a file to write interval records to
        :type filename: str
        """
        self.__filename = filename
        self.__seq = ''

    def __enter__(self):
        self.__output = open(self.__filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()

    def write(self, interval_record):
        """
        Given an interval format record, write it to the file
        specified when the object was created.

        :param interval_record: an interval format record to be
            written to the file
        :type interval_record: Record
        """
        if self.__seq != interval_record.seq:
            self.__seq = interval_record.seq
            self.__output.write('>{}\n'.format(self.__seq))
        self.__output.write('{} - {}\n'.format(interval_record.start,
                                               interval_record.end))
