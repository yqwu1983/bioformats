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
