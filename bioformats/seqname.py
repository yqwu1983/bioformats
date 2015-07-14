#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
from future.utils import iteritems
from .exception import SeqRenameError
from .exception import IncorrectDictError
from .exception import MissingSeqNameError

logging.basicConfig()
logger = logging.getLogger(__name__)


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
            lineno = 0
            for line in dict_file:
                line = line.rstrip()
                lineno += 1
                if line:
                    words = line.rstrip().split(None, 2)
                    if len(words) < 2:
                        logger.error('line %d: incorrect renaming '
                                     'dictionary entry %s', lineno,
                                     line)
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

    def renamed(self, filename, reverse=False, no_desc=False):
        """
        The method implements a generator of lines from the specified
        FASTA file with changed sequence names.

        :param filename: a name of a FASTA file which sequences are 
        to be renamed
        :param reverse: if specified, new sequence names will be 
        replaced with the original ones
        :param no_desc: remove a FASTA sequence description
        :type filename: str
        :type reverse: bool
        :type no_desc: bool
        :return: a line from a FASTA file with renamed sequences
        :rtype: str
        """
        # select the renaming dictionary
        renaming_dict = self.reverse_renaming_dict if reverse else \
            self.renaming_dict
        with open(filename) as fasta_file:
            lineno = 0
            for line in fasta_file:
                lineno += 1
                if line.startswith('>'):
                    # the line is a sequence header, rename it
                    line_parts = line.split()
                    seq_name = line_parts[0][1:]
                    if seq_name not in renaming_dict:
                        logger.error('line %d: missing sequence %s',
                                     lineno, seq_name)
                        raise MissingSeqNameError(seq_name)
                    else:
                        line_parts[0] = '>' + renaming_dict[seq_name]
                    if no_desc:
                        yield line_parts[0] + '\n'
                    else:
                        yield ' '.join(line_parts) + '\n'
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
            lineno = 0
            for line in table_file:
                line = line.rstrip()
                lineno += 1
                if line.startswith(comment_char):
                    yield line + '\n'
                elif line:
                    line_parts = line.split(sep)
                    if line_parts[column] not in renaming_dict:
                        logger.error('line %d: missing sequence %s',
                                     lineno, line_parts[column])
                        raise MissingSeqNameError(
                            line_parts[column])
                    else:
                        line_parts[column] = renaming_dict[
                            line_parts[column]]
                    yield sep.join(line_parts) + '\n'


class NcbiBaseSeqRenamer(BaseSeqRenamer):
    """
    This abstract class adds routines to process NCBI accession number
    files to the BaseSeqRenamer class.
    """

    acceptable_formats = ('chr',
                          'refseq_full', 'genbank_full',
                          'refseq_gi', 'genbank_gi',
                          'refseq', 'genbank',
                          'chr_refseq', 'chr_genbank')

    @staticmethod
    def form_identifier(fmt, chrom, refseq, refseq_gi, genbank,
                        genbank_gi, remove_seq_version=False):
        """
        Given a format and data from an NCBI accession number file,
        return a sequence name in the specified format.

        :param fmt: a format of the required sequence names, possible
            values: chr, refseq_full, genbank_full, refseq_gi,
            genbank_gi, refseq_acc_num, genbank_acc_num
        :param chrom: a chromosome or fragment name
        :param refseq: a RefSeq accession number
        :param refseq_gi: a gene identifier in RefSeq
        :param genbank: a GenBank accession number
        :param genbank_gi: a gene identifier in GenBank
        :param remove_seq_version: remove sequence version from an
            accession number
        :type fmt: str
        :type chrom: str
        :type refseq: str
        :type refseq_gi: str
        :type genbank: str
        :type genbank_gi: str
        :type remove_seq_version: bool
        :return: a sequence name in the specified format
        :rtype: str
        """
        # if required, remove sequence version from accession numbers
        if remove_seq_version:
            refseq = refseq.split('.')[0]
            genbank = genbank.split('.')[0]

        if fmt == 'refseq_full':
            result = 'gi|' + refseq_gi + '|ref|' + refseq + '|'
        elif fmt == 'genbank_full':
            result = 'gi|' + genbank_gi + '|gb|' + genbank + '|'
        elif fmt == 'refseq_gi':
            result = refseq_gi
        elif fmt == 'genbank_gi':
            result = genbank_gi
        elif fmt == 'refseq':
            result = refseq
        elif fmt == 'genbank':
            result = genbank
        elif fmt == 'chr_genbank':
            result = chrom + '_' + genbank
        elif fmt == 'chr_refseq':
            result = chrom + '_' + refseq
        else:
            result = chrom

        return result

    def read_ncbi_acc_num(self, filename, orig_fmt, new_fmt, prefix='',
                          suffix='', remove_seq_version=False):
        """
        Given a name of a file with accession numbers obtained from
        NCBI, form a renaming dictionary of specified original and
        new sequence names.

        :param filename: a name of an NCBI file with accession numbers
        :param orig_fmt: a specification of the format of original
            sequence names
        :param new_fmt: a specification of the format of new sequence
            names
        :param prefix: a prefix to be added to new sequence names
        :param suffix: a suffix to be added to new sequence names
        :param remove_seq_version: remove a sequence version from an
            accession number
        :type filename: str
        :type orig_fmt: str
        :type new_fmt: str
        :type prefix: str
        :type suffix: str
        :type remove_seq_version: bool
        """
        # check if correct format values are specified
        if orig_fmt not in NcbiBaseSeqRenamer.acceptable_formats:
            logger.error('incorrect format %s', orig_fmt)
            raise SeqRenameError(orig_fmt)
        if new_fmt not in NcbiBaseSeqRenamer.acceptable_formats:
            logger.error('incorrect format %s', new_fmt)
            raise SeqRenameError(new_fmt)

        with open(filename) as ncbi_acc_num_file:
            for line in ncbi_acc_num_file:
                line = line.rstrip()
                # skip comments
                if line.startswith('#'):
                    continue
                line_parts = line.split(None, 5)
                if len(line_parts) < 5:
                    logger.error('incorrect dictionary entry %s', line)
                    raise IncorrectDictError(line)
                chrom, refseq, refseq_gi, genbank, genbank_gi = \
                    line_parts
                # form a key and a value of the renaming dictionary
                # according to the specified format arguments
                key = NcbiBaseSeqRenamer.form_identifier(orig_fmt,
                                                         chrom,
                                                         refseq,
                                                         refseq_gi,
                                                         genbank,
                                                         genbank_gi)

                value = NcbiBaseSeqRenamer.form_identifier(
                    new_fmt, chrom, refseq, refseq_gi, genbank,
                    genbank_gi, remove_seq_version)

                value = prefix + value + suffix

                self.renaming_dict[key] = value
                self.reverse_renaming_dict[value] = key


class NcbiFastaSeqRenamer(FastaSeqRenamer, NcbiBaseSeqRenamer):
    """
    The class implements functionality of a sequence renamer for a
    FASTA file that can handle NCBI accession numbers.
    """
    pass


class NcbiTableSeqRenamer(TableSeqRenamer, NcbiBaseSeqRenamer):
    """
    The class implements functionality of a sequence renamer for a
    plain-text file of separated values that can handle NCBI accession
    numbers.
    """
    pass
