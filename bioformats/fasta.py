#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import random
import vcf
import pyfaidx
from builtins import range  # pylint:disable=redefined-builtin
from .exception import BioformatsError
from . import bed


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


class FlankNFilter(object):
    """
    The class implements routines for filtering features from a FASTA
    or VCF file by the presence of N's in their neighborhood.
    """
    def __init__(self, fasta_filename, flank_len=100):
        """
        Create a Neighborhood filter object.

        :param fasta_filename: a name of a FASTA file
        :param flank_len: the length of a region to be checked
        :type fasta_filename: str
        :type flank_len: int
        """
        self.__len = flank_len
        self.__fasta = pyfaidx.Fasta(fasta_filename)

    def check_seq(self, seq_id, start_pos, end_pos):
        """
        Check if the specified region contains N's. Positions are
        zero-based and half-opened as in the BED format,

        :param seq_id: a sequence ID
        :param start_pos: a start position in bp
        :param end_pos: a start position in bp
        :type seq_id: str
        :type start_pos: int
        :type end_pos: int
        :return: True if the specified region contains no N's,
            False otherwise
        :rtype: bool
        """
        seq = str(self.__fasta[seq_id][start_pos:end_pos])
        return 'N' not in seq

    def filter_bed(self, input_filename, output_filename, strict):
        """
        Filter routines from the specified BED file and output the
        filtered ones to the specified output file.

        :param input_filename: a name of a BED file containing
            features to be filtered
        :param output_filename: the output BED file name
        :param strict: if specified, require the flank regions to
            have exactly the specified length
        :type input_filename: str
        :type output_filename: str
        :type strict: bool
        """
        with open(input_filename) as input_file:
            bed_reader = bed.Reader(input_file)
            with bed.Writer(output_filename) as output_file:
                for record in bed_reader.records():
                    seq = record.seq
                    seq_len = len(self.__fasta[seq])
                    if strict and record.start < self.__len:
                        continue
                    down_start = max(0, record.start - self.__len)
                    down_end = record.start
                    up_start = record.end
                    if strict and record.end + self.__len > seq_len:
                        continue
                    up_end = min(seq_len, up_start + self.__len)
                    if self.check_seq(seq, down_start, down_end) and \
                            self.check_seq(seq, up_start, up_end):
                        output_file.write(record)

    def filter_vcf(self, input_filename, output_filename, strict):
        """
        Filter routines from the specified BED file and output the
        filtered ones to the specified output file.

        :param input_filename: a name of a BED file containing
            features to be filtered
        :param output_filename: the output BED file name
        :param strict: if specified, require the flank regions to
            have extactly the specified length
        :type input_filename: str
        :type output_filename: str
        :type strict: bool
        """
        with open(input_filename) as input_file:
            vcf_reader = vcf.Reader(input_file)
            with open(output_filename, 'w') as output_file:
                vcf_writer = vcf.Writer(output_file, vcf_reader)
                for record in vcf_reader:
                    seq = record.CHROM
                    seq_len = len(self.__fasta[seq])
                    pos = record.POS - 1
                    if strict and pos < self.__len:
                        continue
                    down_start = max(0, pos - self.__len)
                    down_end = pos
                    up_start = pos + 1
                    if strict and pos + 1 + self.__len > seq_len:
                        continue
                    up_end = min(seq_len, up_start + self.__len)
                    if self.check_seq(seq, down_start, down_end) and \
                            self.check_seq(seq, up_start, up_end):
                        vcf_writer.write_record(record)
