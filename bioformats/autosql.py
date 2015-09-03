#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import re
from collections import namedtuple

logging.basicConfig()
logger = logging.getLogger(__name__)

Table = namedtuple('Table', ('name', 'desc', 'entries'))
TableEntry = namedtuple('TableEntry', ('type', 'num', 'name', 'desc'))


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
            entry_type, entry_name, entry_desc = line.rstrip().split(
                None, 3)
            # check if the field type is an array
            array_matches = re.search('[\d+]', entry_type, 1)
            if array_matches is not None:
                entry_num = int(array_matches.group(0))
            else:
                entry_num = None
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
