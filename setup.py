#!/usr/bin/env python
from setuptools import setup
import sys, os

version = '0.0.1'

setup(name='sufam',
      version=version,
      description="So U Found A Mutation?",
      long_description="""Validate mutations from supplied vcf in supplied bam""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python validation mutation mpileup samtools pileup vcf bam',
      author='Ino de Bruijn',
      author_email='ino@ino.pm',
      maintainer='Ino de Bruijn',
      maintainer_email='ino@ino.pm',
      url='https://github.com/inodb/sufam',
      license='FreeBSD',
      packages=['sufam'],
      include_package_data=True,
      zip_safe=False,
      install_requires=['numpy>=1.9.2',
                        'pandas>=0.16.1',
                        ],
      entry_points={
          'console_scripts': [
              'sufam = sufam.__main__:main'
          ]
      },
      )
