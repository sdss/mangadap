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

_IDLUTILS_VER = 'v5_5_34'
_MANGADRP_VER = 'v2_7_1'
_MANGACORE_VER = 'v1_8_0'

_VERSION = '2.5.2'
_RELEASE = 'dev' not in _VERSION
_MINIMUM_PYTHON_VERSION = '3.5'

def get_package_data(root='python'):
    """Generate the list of package data."""
    return [os.path.relpath(f, root) 
                    for f in glob.glob(os.path.join(root, 'mangadap/config/*/*.ini'))] \
            + [os.path.relpath(f, root) 
                    for f in glob.glob(os.path.join(root, 'mangadap/tests/data/*'))]

def get_data_files(root=['data'], ext=['par', 'fits', 'fits.gz'], depth=2):
    """Generate the list of data files."""
    data_files = []
    for r in root:
        _data_files = []
        for d in range(depth):
            _d = '/'.join(['*']*(d+1))
            for e in ext:
                _data_files += glob.glob(os.path.join(r, _d, '*.{0}'.format(e)))
        data_files += [(r, _data_files)]
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


def run_setup(package_data, data_files, scripts, packages, install_requires):

    setup(name='sdss-mangadap',
          version=_VERSION,
          license='BSD3',
          description='MaNGA Data Analysis Pipeline',
          long_description=open('README.md').read(),
          long_description_content_type='text/markdown',
          author='SDSS-IV/MaNGA Pipeline Group',
          author_email='westfall@ucolick.org',
          keywords='astronomy analysis-pipeline spectroscopy MaNGA',
          url='https://github.com/sdss/mangadap',
          python_requires='>='+_MINIMUM_PYTHON_VERSION,
          packages=packages,
          package_dir={'': 'python'},
          package_data={'': package_data},
          include_package_data=True,
          data_files=data_files,
          install_requires=install_requires,
          scripts=scripts,
          setup_requires=[ 'pytest-runner' ],
          tests_require=[ 'pytest' ],
          classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
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
        if k not in os.environ and 'MANGACORE' not in k:
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

    # Get the package data (data inside the main product root)
    package_data = get_package_data()

    # Compile the additional data files to include (data outside the
    # main product root)
    data_files = get_data_files()

    # Compile the scripts in the bin/ directory
    scripts = get_scripts()

    # Get the packages to include
    packages = find_packages(where='python')

    # Collate the dependencies based on the system text file
    install_requires = get_requirements()

    # Run setup from setuptools
    run_setup(package_data, data_files, scripts, packages, install_requires)

    # Check if the environmental variables are found and warn the user
    # of their defaults
    check_environment()

