#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
from . import bed
from collections import namedtuple
from .exception import RepeatMaskerError

logging.basicConfig()
logger = logging.getLogger(__name__)

repeatmasker_out_columns = ('sw_score', 'subst_perc', 'del_perc',
                            'ins_perc', 'query', 'query_start',
                            'query_end', 'query_past',
                            'is_complement', 'repeat_name',
                            'repeat_class_family', 'repeat_prior',
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

        # process the last optional column
        if line_parts[-1].endswith('*'):
            line_parts[-1] = line_parts[-1].split(' ', 1)[0]
            line_parts += ['*']
        else:
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


repeat_class_colors = {
    'LINE': (0, 0, 0),      # black
    'SINE': (255, 0, 0),    # red
    'LTR': (0, 205, 0),     # green3
    'DNA': (0, 0, 255),     # blue
    'Simple_repeat': (0, 255, 255),     # cyan
    'Low_complexity': (255, 0, 255),    # magenta
    'Satellite': (255, 255, 0),     # yellow
    'RNA': (190, 190, 190),         # gray
    'Other': (30, 144, 255),        # dodgerblue
    'Unknown': (178, 34, 34)        # firebrick
    }


def whiten_color(color, alpha):
    """
    Given an RGB color, whiten it to the specified alpha.

    :param color: an RGB color as a tuple of three numbers
    :param alpha: a whitening alpha from 0 to 100
    :type color: tuple
    :type alpha: float
    :return: a whitened color
    :rtype: tuple
    """
    if not (0 <= alpha <= 100):
        logger.error('incorrect color whitening alpha %d', alpha)
        raise RepeatMaskerError
    return tuple(int(round(x + (255 - x)/100.0*alpha)) for x in color)


def rmout2bed_record(rm_record, name='id', color='class',
                     is_short=False):
    """
    Convert a record from a RepeatMasker out file to a BED record.

    :param rm_record: a record from a RepeatMasker out file
    :param name: how to form a BED record name
    :param color: how to form a BED record color
    :param is_short: if True, then only repeat coordinates are given
        in the BED record
    :type rm_record: Record
    :param name: str
    :param color: str
    :param is_short: bool
    :return: the BED record corresponding to the specified
        RepeatMasker out record
    :rtype: bed.Record
    """
    if is_short:
        bed_record = bed.Record(
            seq=rm_record.query,
            start=rm_record.query_start - 1,
            end=rm_record.query_end,
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
    else:
        # choose a BED record name
        if name == 'id':
            bed_name = 'ID_' + rm_record.id
        elif name == 'name':
            bed_name = rm_record.repeat_name
        elif name == 'class':
            bed_name = rm_record.repeat_class_family.split('/')[0]
        elif name == 'family':
            bed_name = rm_record.repeat_class_family.split('/')[1]
        elif name == 'class_family':
            bed_name = rm_record.repeat_class_family
        else:
            logger.error('incorrect name parameter %s', name)
            raise RepeatMaskerError

        # choose a BED record color
        repeat_class = rm_record.repeat_class_family.split('/')[0]
        # remove question marks if there are any
        repeat_class = repeat_class.strip('?')
        # process various RNA classes
        if 'RNA' in repeat_class and repeat_class != 'RNA':
            repeat_class = 'RNA'
        repeat_identity = int(round(10 * (100 - float(
            rm_record.subst_perc))))
        if repeat_class not in repeat_class_colors:
            logger.info('classifying repeat class %s as unknown',
                        repeat_class)
            repeat_class = 'Unknown'
        if color == 'class':
            bed_color = repeat_class_colors[repeat_class]
        elif color == 'identity':
            bed_color = whiten_color((0, 0, 0), 100 -
                                     repeat_identity/10)
        elif color == 'class_identity':
            bed_color = whiten_color(repeat_class_colors[repeat_class],
                                     100 - repeat_identity/10)
        else:
            logger.error('incorrect color parameter %s', color)
            raise RepeatMaskerError

        bed_record = bed.Record(
            seq=rm_record.query,
            start=rm_record.query_start - 1,
            end=rm_record.query_end,
            name=bed_name,
            score=repeat_identity,
            strand='-' if rm_record.is_complement == 'C' else '+',
            thick_start=rm_record.query_start,
            thick_end=rm_record.query_end,
            color=','.join(map(str, bed_color)),
            block_num=None,
            block_sizes=None,
            block_starts=None,
            extra=[
                rm_record.sw_score,
                rm_record.subst_perc,
                rm_record.del_perc,
                rm_record.ins_perc,
                rm_record.query_past,
                rm_record.repeat_name,
                rm_record.repeat_class_family,
                rm_record.repeat_prior,
                rm_record.repeat_start,
                rm_record.repeat_end,
                rm_record.id
            ]
        )
    return bed_record
