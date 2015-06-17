#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from setuptools import setup, find_packages

setup(name='bioformats',
      version='0.1.0.dev1',
      description='Classes to handle bioinformatics data',

      author='Gaik Tamazian',
      author_email='gaik.tamazian@gmail.com',

      license='MIT',

      classifiers=[
          "Development Status :: 2 - Pre-Alpha",
          "License :: OSI Approved :: MIT License",
          "Environment :: Console",
          "Interned Audience :: Science/Research",
          "Natural Language :: English",
          "Operating System :: Unix",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.4",
          "Topic :: Scientific/Engineering :: Bio-Informatics"
      ],

      packages=find_packages(exclude=['doc', 'tests*']),

      install_requires=['pyfaidx',
                        'future',
                        'six'],
      )
