# !usr/bin/env python
# -*- coding: utf-8 -*-
# K. Westfall, 12 May 2018
#   Adapted from Marvin's setup.py file

import sys
import os
import glob

from setuptools import setup, find_packages

from astropy.utils.data import download_file

import argparse
import shutil


def get_data_files():
    """Generate the list of data files."""

    data_files = []
    data_roots = [ 'data', 'python/mangadap/config' ]
    for root in data_roots:
        for path, directories, files in os.walk(root):
            for f in files:
                data_path = '/'.join(path.split('/')[1:])
                data_files.append(os.path.join(data_path, f))
    return data_files


def get_scripts():
    """ Grab all the scripts in the bin directory.  """
    scripts = []
    if os.path.isdir('bin'):
        scripts = [ fname for fname in glob.glob(os.path.join('bin', '*'))
                                if not os.path.basename(fname).endswith('.rst') ]
    return scripts


def get_requirements():
    """ Get the package requirements from the system file. """
    requirements_file = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
    return install_requires



# requirements
requirements_file = os.path.join(os.path.dirname(__file__), 'requirements.txt')
install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                    if not line.strip().startswith('#') and line.strip() != '']


NAME = 'sdss-mangadap'
# do not use x.x.x-dev.  things complain.  instead use x.x.xdev
VERSION = '2.2.2dev'
RELEASE = 'dev' not in VERSION
minimum_python_version = '3.5'

def run_setup(data_files, scripts, packages, install_requires):

    setup(name=NAME,
          version=VERSION,
          license='BSD3',
          description='MaNGA Data Analysis Pipeline',
          long_description=open('README.md').read(),
          author='SDSS-IV/MaNGA Pipeline Group',
          author_email='westfall@ucolick.org',
          keywords='astronomy MaNGA',
          url='https://github.com/sdss/mangadap',
          python_requires='>='+minimum_python_version,
          packages=packages,
          package_dir={'': 'python'},
          package_data={'': data_files},
          include_package_data=True,
          install_requires=install_requires,
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3 :: Only'
              'Topic :: Documentation :: Sphinx',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Topic :: Software Development :: Libraries :: Python Modules'
          ],
          )

#-----------------------------------------------------------------------
if __name__ == '__main__':

    # Compile the data files to include
    data_files = get_data_files()
    # Compile the scripts in the bin/ directory
    scripts = get_scripts()
    # Get the packages to include
    packages = find_packages()
    # Collate the dependencies based on the system text file
    install_requires = get_requirements()
    # Run setup from setuptools
    run_setup(data_files, scripts, packages, install_requires)

