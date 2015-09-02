#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import pyfaidx
import random
from builtins import range  # pylint:disable=redefined-builtin
from .exception import BioformatsError


class Writer(object):
    """
    The class implements routines to write sequences in the FASTA
    format.
    """

    def __init__(self, filename, width=72):
        """
        Create a Writer object to write sequences in a FASTA file.

        :param filename: a name of a file to write sequences to
        :type filename: str
        """
        self.__filename = filename
        self.__width = width

    def __enter__(self):
        self.__output = open(self.__filename, 'w')
        return self

    def write(self, header, sequence):
        """
        Write a sequence to the FASTA file.

        :param header: a sequence header
        :param sequence: a sequence
        :type header: str
        :type sequence: str
        """
        self.__output.write('>{}\n'.format(header))
        seq_lines = []
        for i in range(0, len(sequence), self.__width):
            seq_lines.append(sequence[i:i + self.__width])
        self.__output.write('\n'.join(seq_lines))
        self.__output.write('\n')

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()


class RandomSequence(object):
    """
    The class implements routines to create random nucleotide
    sequences.
    """

    def __init__(self, length):
        """
        Create an object to generate random nucleotide sequences of
        the specified length.

        :param length: a sequence length
        :type length: int
        """
        self.__length = length

    def get(self):
        """
        Get a random nucleotide sequence.

        :return: a random nucleotide sequence
        :rtype: str
        """
        result = []
        nucleotides = ('A', 'C', 'G', 'T')
        for _ in range(self.__length):
            result.append(random.choice(nucleotides))

        return ''.join(result)


class Reorder(object):
    """
    The class implements reordering of sequences in a FASTA file.
    """

    def __init__(self, order_filename):
        """
        Given a name of a file containing sequence names in the
        specific order, read this order.

        :param order_filename: a name of a file with ordered sequence
            names
        """
        self.__order = []
        with open(order_filename) as order_file:
            for line in order_file:
                self.__order.append(line.rstrip())

    @property
    def order(self):
        return self.__order

    def write(self, input_file, output_file, ignore_missing=True):
        """
        Given a handle of an input FASTA file, reorder its sequences
        and write to the specified output file.

        :param input_file: a name of an input FASTA file
        :param output_file: a name of an output FASTA file
        :param ignore_missing: ignore sequences given in the sequence
            order but present in the input FASTA file
        :type ignore_missing: bool
        """
        seq_reader = pyfaidx.Fasta(input_file)
        with Writer(output_file) as seq_writer:
            for i in self.__order:
                if i not in seq_reader.keys():
                    if not ignore_missing:
                        raise BioformatsError('missing sequence {'
                                              '}'.format(i))
                else:
                    seq_writer.write(i, str(seq_reader[i]))
