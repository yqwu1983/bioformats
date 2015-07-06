#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from future.utils import iteritems
from .exception import IncorrectDictError
from .exception import MissingSeqNameError


class BaseSeqRenamer(object):
    """
    This class implements basic functionality of a sequence renamer
    class: storing renaming dictionaries (direct and reverse),
    reading and writing them to and from a file.
    """
    def __init__(self):
        """
        Create a BaseSeqRenamer object.
        """
        self.renaming_dict = dict()
        self.reverse_renaming_dict = dict()

    def read_renaming_dict(self, filename):
        """
        Given a name of a file, read the renaming dictionary from it.
        Note that new entries are added to the existing ones.

        :param filename: a name of a file containing a renaming
            dictionary
        :type filename: str
        """
        with open(filename) as dict_file:
            for line in dict_file:
                line = line.rstrip()
                if line:
                    words = line.rstrip().split(None, 2)
                    if len(words) < 2:
                        raise IncorrectDictError(line)
                    orig_name, new_name = words[0], words[1]
                    self.renaming_dict[orig_name] = new_name
                    self.reverse_renaming_dict[new_name] = orig_name

    def write_renaming_dict(self, filename):
        """
        Given a name of a file, write the renaming dictionary to it.

        :param filename: a name of a file to write the renaming
            dictionary to
        :type filename: str
        """
        with open(filename, 'w') as dict_file:
            for orig_name, new_name in iteritems(self.renaming_dict):
                dict_file.write(orig_name + '\t' + new_name + '\n')


class FastaSeqRenamer(BaseSeqRenamer):
    """
    This class implements functionality of a sequence renamer for a
    FASTA file. Names of sequences in a FASTA file are replaced with
    the corresponding values from the renaming dictionary.
    """

    def renamed(self, filename, reverse=False):
        """
        The method implements a generator of lines from the specified
        FASTA file with changed sequence names.

        :param filename: a name of a FASTA file which sequences are 
        to be renamed
        :param reverse: if specified, new sequence names will be 
        replaced with the original ones
        :type filename: str
        :type reverse: bool
        :return: a line from a FASTA file with renamed sequences
        :rtype: str
        """
        # select the renaming dictionary
        renaming_dict = self.reverse_renaming_dict if reverse else \
            self.renaming_dict
        with open(filename) as fasta_file:
            for line in fasta_file:
                line = line.rstrip()
                if line.startswith('>'):
                    # the line is a sequence header, rename it
                    line_parts = line.split()
                    seq_name = line_parts[0][1:]
                    if seq_name not in renaming_dict:
                        raise MissingSeqNameError(seq_name)
                    else:
                        line_parts[0] = '>' + renaming_dict[seq_name]
                    yield ' '.join(line_parts)
                else:
                    # the line contains a sequence, just return it
                    yield line


class TableSeqRenamer(BaseSeqRenamer):
    """
    This class implements functionality of a sequence renamer for a
    plain-text file of separated values.
    """

    def renamed(self, filename, column, reverse=False, sep='\t',
                comment_char='#'):
        """
        The method implements a generator of lines from the specified
        plain-text file of separated values.

        :param filename: a name of a file which lines are to be
        searched for sequence names to be changed
        :param column: the number of a column which contains sequence
        names to be changed
        :param reverse: if specified, new sequence names will be
        replaced with the original ones
        :param sep: a separator symbol which delimits columns in the
        specified file
        :param comment_char: a character designating a comment line
        :type filename: str
        :type column: int
        :type reverse: bool
        :type sep: str
        :type comment_char: str
        :return: a line from the specified file with renamed sequences
        in the specified columns
        :rtype: str
        """
        # select the renaming dictionary
        renaming_dict = self.reverse_renaming_dict if reverse else \
            self.renaming_dict
        with open(filename) as table_file:
            for line in table_file:
                line = line.rstrip()
                if line.startswith(comment_char):
                    yield line
                elif line:
                    line_parts = line.split(sep)
                    if line_parts[column] not in renaming_dict:
                        raise MissingSeqNameError(
                            line_parts[column])
                    else:
                        line_parts[column] = renaming_dict[
                            line_parts[column]]
                    yield sep.join(line_parts)
