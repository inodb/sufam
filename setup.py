#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
import versioneer

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

setup(name='sufam',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="So U Found A Mutation?",
      long_description=read("README.rst"),
      classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      keywords='Python validation mutation mpileup samtools pileup vcf bam',
      author='Ino de Bruijn',
      author_email='ino@ino.pm',
      url='https://github.com/inodb/sufam',
      license="MIT",
      include_package_data=True,
      packages=find_packages(exclude=['test*']),
      install_requires=['numpy>=1.9.2',
                        'pandas>=0.16.1',
                        ],
      entry_points={
          'console_scripts': [
              'sufam = sufam.__main__:main'
          ]
      },
      )
