#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple
from .exception import IntervalError
from . import bed
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
    def __init__(self, handle):
        """
        Given a handle of a file, create an interval format reader
        object to read data from it.

        :param handle: a handle of an interval format file
        """
        self.__handle = handle
        self.__line = ''
        self.__lineno = 0
        self.__seq = ''

    def intervals(self):
        """
        Iterate through records in the interval format file the
        object was created from

        :return: a record from the interval format file the object
            was created from
        :rtype: Record
        """
        for self.__line in self.__handle:
            self.__line = self.__line.rstrip()
            self.__lineno += 1
            if self.__line.startswith('>'):
                # remove the comment if present
                self.__line = self.__line.split()[0]
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


def convert_interval2bed(interval_filename, bed_filename):
    """
    Convert a specified interval file to the BED format.

    :param interval_filename: a name of an interval file
    :param bed_filename: a name of the output BED file
    :type interval_filename: str
    :type bed_filename: str
    """
    with open(interval_filename) as input_file:
        reader = Reader(input_file)
        with bed.Writer(bed_filename) as writer:
            for record in reader.intervals():
                bed_record = bed.Record(
                    seq=record.seq,
                    start=record.start - 1,
                    end=record.end,
                    name=None,
                    score=None,
                    strand=None,
                    thick_start=None,
                    thick_end=None,
                    color=None,
                    block_num=None,
                    block_sizes=None,
                    block_starts=None,
                    extra=[]
                )
                writer.write(bed_record)
