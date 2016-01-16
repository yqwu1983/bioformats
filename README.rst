|Travis| |PyPI| |Landscape| |Coveralls|

==========
bioformats
==========

A collection of routines and related Python classes to handle
bioinformatics data.

Installation
------------

We recommend to install bioformats via **pip**::

    pip install bioformats

**pip** automatically resolves bioformats dependencies and installs
required packages.

Usage
-----

The package tool are called via the main wrapper script **bioformats**.
For example, one may get the genotype count table from a VCF file by
calling the **vcfgeno2bed** tool in the following way.

    bioformats vcfgeno2bed variants.vcf genotype_counts.txt

For every tool, the help message is available by using the **-h**
option. Also one may get the list of the package tools by calling
**bioformats** with the **-h** option.

.. |PyPI| image:: https://img.shields.io/pypi/v/bioformats.svg?branch=master
    :target: https://pypi.python.org/pypi/bioformats
.. |Travis| image:: https://travis-ci.org/gtamazian/bioformats.svg?branch=master
    :target: https://travis-ci.org/gtamazian/bioformats
.. |Coveralls| image:: https://coveralls.io/repos/gtamazian/bioformats/badge.svg?branch=master 
    :target: https://coveralls.io/r/gtamazian/bioformats?branch=master
.. |Landscape| image:: https://landscape.io/github/gtamazian/bioformats/master/landscape.svg?style=flat
   :target: https://landscape.io/github/gtamazian/bioformats/master
   :alt: Code Health
