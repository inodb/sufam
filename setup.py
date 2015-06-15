#!/usr/bin/env python
from setuptools import setup
import sys, os
import versioneer

setup(name='sufam',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="So U Found A Mutation?",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python validation mutation mpileup samtools pileup vcf bam',
      author='Ino de Bruijn',
      author_email='ino@ino.pm',
      maintainer='Ino de Bruijn',
      maintainer_email='ino@ino.pm',
      url='https://github.com/inodb/sufam',
      download_url = 'https://github.com/inodb/sufam/tarball/'+versioneer.get_version(),
      packages=['sufam'],
      install_requires=['numpy>=1.9.2',
                        'pandas>=0.16.1',
                        ],
      entry_points={
          'console_scripts': [
              'sufam = sufam.__main__:main'
          ]
      },
      )
