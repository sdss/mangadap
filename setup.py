#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# K. Westfall, 22 May 2018
#   Adapted from Marvin's setup.py file

# Imports
import os
import glob
from setuptools import setup, find_packages

import requests
import warnings

_IDLUTILS_VER = 'v5_5_30'
_MANGADRP_VER = 'v2_4_3'
_MANGACORE_VER = 'v1_6_2'

_VERSION = '2.2.3dev'
_RELEASE = 'dev' not in _VERSION
_MINIMUM_PYTHON_VERSION = '3.5'

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
                                if not os.path.basename(fname).endswith('.rst') and
                                   not os.path.basename(fname).endswith('.bash') ]
    return scripts


def get_requirements():
    """ Get the package requirements from the system file. """
    requirements_file = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
#    install_requires = [line.strip() for line in open(requirements_file)
#                        if not line.strip().startswith('#') and line.strip() != '']
    return install_requires


def run_setup(data_files, scripts, packages, install_requires):

    setup(name='sdss-mangadap',
          version=_VERSION,
          license='BSD3',
          description='MaNGA Data Analysis Pipeline',
          long_description=open('README.md').read(),
          author='SDSS-IV/MaNGA Pipeline Group',
          author_email='westfall@ucolick.org',
          keywords='astronomy analysis-pipeline spectroscopy MaNGA',
          url='https://github.com/sdss/mangadap',
          python_requires='>='+_MINIMUM_PYTHON_VERSION,
          packages=packages,
          package_dir={'': 'python'},
          package_data={'': data_files},
          include_package_data=True,
          install_requires=install_requires,
          scripts=scripts,
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3.5',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3 :: Only',
              'Topic :: Documentation :: Sphinx',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Topic :: Software Development :: Libraries :: Python Modules'
          ])

# TODO: Put these in a config file
def default_paths():
    return { 'MANGADRP_VER': _MANGADRP_VER,
             'MANGA_SPECTRO_REDUX': os.path.join(os.environ['HOME'], 'MaNGA', 'redux'),
             'MANGADAP_VER': _VERSION,
             'MANGA_SPECTRO_ANALYSIS': os.path.join(os.environ['HOME'], 'MaNGA', 'analysis'),
             'MANGACORE_VER': _MANGACORE_VER,
             'MANGACORE_DIR': os.path.join(os.environ['HOME'], 'MaNGA', 'core', _MANGACORE_VER)
           }


def check_environment():
    ev = default_paths()
    for k in ev.keys():
        if k not in os.environ:
            warnings.warn('{0} environmental variable undefined.  Default is: {1}'.format(k,ev[k]))


def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s\n' % (category.__name__, message)


#-----------------------------------------------------------------------
if __name__ == '__main__':

    warnings.formatwarning = short_warning

    # Pull over the maskbits file
    idlpath = 'https://svn.sdss.org/public/repo/sdss/idlutils'
    ofile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'sdss',
                         'sdssMaskbits.par')
    try:
        idldir = 'tags/{0}'.format(_IDLUTILS_VER)
        parpath = '{0}/{1}/data/sdss/sdssMaskbits.par'.format(idlpath,idldir)
        r = requests.get(parpath)
        open(ofile, 'wb').write(r.content)
    except:
        try:
            idldir = 'trunk'
            parpath = '{0}/{1}/data/sdss/sdssMaskbits.par'.format(idlpath,idldir)
            open(ofile, 'wb').write(requests.get(parpath).content)
        except:
            warnings.warn('Could not download SDSS maskbits file!')

    # Compile the data files to include
    data_files = get_data_files()

    # Compile the scripts in the bin/ directory
    scripts = get_scripts()

    # Get the packages to include
    packages = find_packages(where='python')

    # Collate the dependencies based on the system text file
    install_requires = get_requirements()

    # Run setup from setuptools
    run_setup(data_files, scripts, packages, install_requires)

    # Check if the environmental variables are found and warn the user
    # of their defaults
    check_environment()

