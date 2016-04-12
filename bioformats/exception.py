#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


class BioformatsError(Exception):
    """
    The class describes a basic error that may occur in routines of
    the bioformats package.
    """
    pass


class BedError(BioformatsError):
    """
    The class describes an error related to the BED format parser.
    """
    pass


class BlastTabError(BioformatsError):
    """
    The class describes an error related to the BLAST tabular format
    parser.
    """
    pass


class LavError(BioformatsError):
    """
    The class describes an error related to the LAV alignment format
    parser.
    """
    pass


class Gff3Error(BioformatsError):
    """
    The class describes an error related to the GFF3 format parser.
    """
    pass


class VcftoolsError(BioformatsError):
    """
    This class describes an error related to any parser of VCFtools
    output files.
    """
    pass


class FrqCountReaderError(VcftoolsError):
    """
    This class describes an error related to the frequency count
    format of VCFtools.
    """
    pass


class IntervalError(BioformatsError):
    """
    This class describes an error related to the interval format used
    by WindowMasker and DustMasker.
    """
    pass


class SeqRenameError(BioformatsError):
    """
    This class describes an error related to sequence renaming
    routines.
    """
    pass


class IncorrectDictError(SeqRenameError):
    """
    This class corresponds to a situation when a specified renaming
    dictionary file contains errors, for example, a column of new names
    is missing.
    """
    pass


class MissingSeqNameError(SeqRenameError):
    """
    This class describes an exception to be raised when a sequence name
    is missing from the dictionary of a sequence renamer object.
    """
    pass


class RepeatMaskerError(BioformatsError):
    """
    This class describes an error related to the RepeatMasker out
    file parser.
    """
    pass


class SnpEffError(BioformatsError):
    """
    This class describes an error related to the snpEff annotation
    parsing.
    """
    pass


class AgpError(BioformatsError):
    """
    This class describes an error related to the AGP processing
    routines.
    """
    pass
