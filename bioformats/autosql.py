#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import re
from collections import namedtuple
from future.utils import iteritems

logging.basicConfig()
logger = logging.getLogger(__name__)

Table = namedtuple('Table', ('name', 'desc', 'entries'))
TableEntry = namedtuple('TableEntry', ('type', 'num', 'name', 'desc'))

# below we introduce the hierarchy of autoSql types
type_order = dict([
    ('byte', 1),
    ('ubyte', 2),
    ('short', 3),
    ('ushort', 4),
    ('int', 5),
    ('uint', 6),
    ('float', 7),
    ('string', 8),
    ('lstring', 9)
])

type_names = {value: key for key, value in iteritems(type_order)}

signed_types = ('byte', 'short', 'int')
unsigned_types = ('ubyte', 'ushort', 'uint')


class Reader(object):
    """
    This class provides routines to read autoSql files.
    """
    def __init__(self, handle):
        """
        Create a Reader object to read entries from an autoSql file.
        """
        self.__handle = handle
        # obtain a table name and desciption
        line = self.__handle.readline().rstrip()
        self.__name = line.split()[1]
        line = self.__handle.readline().rstrip()
        self.__desc = line.strip('"')

        # read a line denoting the start of entry records
        self.__handle.readline().rstrip()

        logging.debug('started reading table %s (%s)...', self.__name,
                      self.__desc)

    @property
    def table_name(self):
        return self.__name

    @property
    def table_desc(self):
        return self.__desc

    def entries(self):
        """
        Return an iterator to iterate through entries of an autoSql
        file.

        :return: an itereator to iterate through autoSql file entries
        """
        for line in self.__handle:
            line = line.rstrip()
            if line != ')':
                entry_type, entry_name, entry_desc = line.split(None, 2)
                # check if the field type is an array
                array_matches = re.search('[\d+]', entry_type, 1)
                if array_matches is not None:
                    entry_num = int(array_matches.group(0))
                    entry_type = entry_type[:array_matches.start(0)-1]
                else:
                    entry_num = None
                entry_name = entry_name.rstrip(';')
                entry_desc = entry_desc.strip('"')
                yield TableEntry(entry_type, entry_num, entry_name,
                                 entry_desc)

        logging.debug('finished reading table %s (%s)', self.__name,
                      self.__desc)

    def get_table(self):
        """
        Return the whole table from an autoSql file at once.

        :return: the tuple representing contents of an autoSql file
        :rtype: Table
        """
        entries = []
        for i in self.entries():
            entries.append(i)
        return Table(self.table_name, self.table_desc, entries)


class Writer(object):
    """
    This class provides routines to write autoSql files.
    """

    def __init__(self, filename, table_name, table_desc):
        """
        Create a Writer object to write entries to an autoSql file.

        :param filename: a name of a file to write entries to
        :param table_name: a name of a table
        :param table_desc: a description of a table
        :type filename: str
        """
        self.__filename = filename
        self.__tbl_name = table_name
        self.__tbl_desc = table_desc

    def __enter__(self):
        logger.debug('started writing table %s(%s) to %s...',
                     self.__tbl_name, self.__tbl_desc, self.__filename)

        self.__output = open(self.__filename, 'w')
        self.__output.write('table {}\n"{}"\n'.format(self.__tbl_name,
                                                      self.__tbl_desc))
        self.__output.write('(\n')
        return self

    def write(self, entry):
        """
        Write an entry to a file in the autoSql format.

        :param entry: an entry from an autoSql file
        :type entry: TableEntry
        """
        line = entry.type
        if entry.num is not None:
            line += '[{}]'.format(entry.num)
        line += ' {}; "{}"\n'.format(entry.name, entry.desc)
        self.__output.write(line)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.write(')\n')
        self.__output.close()

        logger.debug('finished writing table %s (%s) to %s',
                     self.__tbl_name, self.__tbl_desc, self.__filename)


def is_int(x):
    """
    Given a string value, determine if it is an integer.

    :param x: a string value to check if it is an integer
    :type x: str
    :return: whether the specified value is an integer or not
    :rtype: bool
    """
    try:
        int(x)
    except ValueError:
        return False
    return True


def get_int_type(x):
    """
    Given an integer value, determine which the least autoSql type
    that suits it.

    :param x: an integer value
    :type x: int
    :return: an autoSql type for the specified value
    :rtype: str
    """
    if x < 0:
        # consider signed types
        if -pow(2, 7) - 1 < x < pow(2, 7):
            return 'byte'
        elif -pow(2, 15) - 1 < x < pow(2, 15):
            return 'short'
        elif -pow(2, 31) - 1 < x < pow(2, 31):
            return 'int'
    else:
        # consider unsigned types
        if x < pow(2, 8):
            if x < pow(2, 7):
                return 'byte'
            else:
                return 'ubyte'
        elif x < pow(2, 16):
            if x < pow(2, 15):
                return 'short'
            else:
                return 'ushort'
        elif x < pow(2, 32):
            if x < pow(2, 31):
                return 'int'
            else:
                return 'uint'
    return None


def is_float(x):
    """
    Given a string value, determine if it is a floating-point number.

    :param x: a string value to check if it is a floating-point number
    :type x: str
    :return: whether the specified value is a floating-point number
        or not
    :rtype: bool
    """
    try:
        float(x)
    except ValueError:
        return False
    return True


def get_autosql_type(value):
    """
    Given a string value, determine the most appropriate autoSql type
    for it.

    :param value: a string value to determine its autoSql type
    :type value: str
    :return: an autoSql type
    :rtype: str
    """
    if is_int(value) and get_int_type(int(value)) is not None:
        value_type = get_int_type(int(value))
    elif is_float(value):
        value_type = 'float'
    elif len(value) < 256:
        value_type = 'string'
    else:
        value_type = 'lstring'

    return value_type


def compare_autosql_types(x, y):
    """
    Given two values of autoSql types, compare them to each other and
    return the type that is more common.

    :param x: the first autoSql type value
    :param y: the second autoSql type value
    :type x: str
    :type y: str
    :return: the most common of the specified autoSql types
    :rtype: str
    """
    if x is None or y is None:
        if x is None:
            return y
        else:
            return x

    is_x_int = x in (signed_types + unsigned_types)
    is_y_int = y in (signed_types + unsigned_types)
    is_x_signed = x in signed_types
    is_y_signed = y in signed_types
    if is_x_int and is_y_int and is_x_signed != is_y_signed:
        # one of the integer types is signed and another one is not
        if is_x_signed:
            signed_arg = type_order[x]
            unsigned_arg = type_order[y]
        else:
            signed_arg = type_order[y]
            unsigned_arg = type_order[x]
        if unsigned_arg < signed_arg:
            result = signed_arg
        else:
            # the largest signed integer type is int which
            # corresponding number is 5
            result = min(max(signed_arg, unsigned_arg - 1) + 2, 5)
    else:
        result = max(type_order[x], type_order[y])

    return type_names[result]


class Classifier(object):
    """
    This class implements routines to assign autoSql types to arrays
    of values.
    """

    def __init__(self):
        """
        Initialize the classifier object.
        """
        self.__data_type = None
        self.__lengths = []

    def add_value(self, value):
        """
        Add a value from the data set which type is being specified.

        :param value: a data set value
        :type value: str
        """
        self.__lengths.append(len(value))
        curr_type = get_autosql_type(value)
        self.__data_type = compare_autosql_types(self.__data_type,
                                                 curr_type)

    def is_array(self):
        """
        Based on lengths of value strings, consider if they can be
        represented as an array.

        :return: if the classified value can be an array of fixed
            length
        :rtype: bool
        """
        return len(set(self.__lengths)) <= 1

    def __check_for_char_array(self):
        """
        Check if the given values form a char array.

        :return: if the given values form a char array
        :rtype: bool
        """
        return self.is_array() and self.__data_type in ('string',
                                                        'lstring')

    @property
    def data_type(self):
        """
        Return the autoSql data type that suits all specified values.

        :return: autoSql data type specifying all values
        :rtype: str
        """
        if self.__check_for_char_array():
            return 'char[{}]'.format(self.__lengths[0])
        else:
            return self.__data_type
