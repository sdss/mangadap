# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class to keep track of MaNGA MPL dependencies and versioning.

.. note::
    This class was created for the survey-level execution of the DAP.
    The dependencies will necessarily grow, and will be handled using
    module files.  So, even as it is, this class has little utility.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/mangampl.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import os
        from distutils.version import StrictVersion

        # For version checking
        import numpy
        import scipy
        import matplotlib
        import astropy
        import pydl

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **2015**: Original Implementation by K. Westfall (KBW)
    | **Apr 2015**: (KBW) Added MPL-3
    | **20 May 2015**: (KBW) Documentation and Sphinx tests. Change
        default to MPL-3
    | **27 Aug 2015**: (KBW) Added (temporary) MPL-4 dependencies
    | **06 Oct 2015**: (KBW) Changed DRP version for MPL-4 to v1_5_0
    | **22 Oct 2015**: (KBW) Changed DRP version for MPL-4 to v1_5_1.  Added
        version dependency for mangacore and sdss_access.
    | **13 May 2016**: (KBW) Added MPL-5 (for now trunk).
    | **10 Oct 2016**: (KBW) Updated for MPL-6 (DAP is trunk
        temporarily).
    | **11 Oct 2016**: (KBW) Moved version checking to here (previously
        in rundap.py).
    | **11 May 2017**: (KBW) Update MPL-6 dependencies; remove
        module_file function

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import os
import warnings
from distutils.version import StrictVersion

# For version checking
import numpy
import scipy
import matplotlib
import astropy
import pydl

# TODO: Add DAP version

class MaNGAMPL:
    """

    Create an object that keeps track of IDLUTILS, SDSS_ACCESS, MANGACORE, and
    MANGADRP dependencies for a given MPL and other intermediate analyses.
    Dependencies are selected by the "MPL" version as provided in the
    following table:

    **SDSS dependencies**
    +--------+-------------+------------+--------+--------+--------+
    |  MPL   | SDSS_ACCESS |  IDLUTILS  |   CORE |   DRP  |   DAP  |
    +========+=============+============+========+========+========+
    | MPL-1  |        None |    v5_5_16 | v1_0_0 | v1_0_0 |   None |
    +--------+-------------+------------+--------+--------+--------+
    | MPL-2  |        None |    v5_5_17 | v1_0_0 | v1_1_2 |   None |
    +--------+-------------+------------+--------+--------+--------+
    | v1_2_0 |        None |    v5_5_19 | v1_1_0 | v1_2_0 |   None |
    +--------+-------------+------------+--------+--------+--------+
    | MPL-3  |        None |    v5_5_22 | v1_1_0 | v1_3_3 | v1_0_0 |
    +--------+-------------+------------+--------+--------+--------+
    | MPL-4  |       0.0.0 |    v5_5_23 | v1_2_0 | v1_5_1 |  1.1.1 |
    +--------+-------------+------------+--------+--------+--------+
    | MPL-5  |       0.2.0 |    v5_5_26 | v1_3_1 | v2_0_1 |  2.0.2 |
    +--------+-------------+------------+--------+--------+--------+
    | MPL-6  |       0.2.1 |    v5_5_29 | v1_5_0 | v2_2_0 |  trunk |
    +--------+-------------+------------+--------+--------+--------+

    **Python versions**
    +--------+--------+--------+--------+------------+---------+-------+----------------+
    | MPL    | python | numpy  |  scipy | matplotlib | astropy |  pydl | Fwd compatible |
    +========+========+========+========+============+=========+=======+================+
    | MPL-5  |  3.5.1 | 1.11.0 | 0.17.1 |       None |   1.1.2 | 0.5.0 |            yes |
    +--------+--------+--------+--------+------------+---------+-------+----------------+
    | trunk  |  3.5.3 | 1.12.0 | 0.19.0 |      2.0.0 |   1.3.0 | 0.5.3 |            N/A |
    +--------+--------+--------+--------+------------+---------+-------+----------------+

    .. note::
        - "Fwd compatible" means that you can run this version of the
          DAP with the most recent version without needed to revert the
          python packages to their previous versions.
        - Only the survey-level routines are dependent on SDSS_ACCESS
          and MANGACORE; the core DAP processing software only depends
          on the listed python packages.

    .. todo::
        - Change python dependencies to "minimum" versions?

    .. warning::

        - There were no dependency of the survey-level DAP code on SDSS_ACCESS
          until MPL-4.
        - For MPL-1, the MANGACORE version was set to trunk.  Here it defaults
          to v1_0_0.
        - Explicit dependence on matplotlib was introduced in MPL-6.

    Args:
        version (str): (**Optional**) Version as selected via matching
            to the "MPL" string in the table listed above. Default:
            MPL-6.
        strictver (bool): (**Optional**) Strictly check the version
            requirements for this version of the dap.  Default is True,
            meaning the code will raise an exception if the expected
            version are not the same as the environment.

    Raises:
        ValueError: Raised if the provided *version* is not one of the
            MPL values in the above table.

    Attributes:
        mplver (str): Version as selected via matching to the "MPL"
            string in the table listed above. Default: MPL-5.
        accessver (str): SDSS_ACCESS version.
        idlver (str): IDLUTILS version.
        corever (str): MANGACORE version.
        drpver (str): MANGADRP version.
        dapver (str): MANGADAP version.
        pythonver (str): Python version
        numpyver (str): Numpy version
        scipyver (str): Scipy version
        matplotlibver (str): Matplotlib version
        astropyver (str): Astropy version
        pydlver (str): Pydl version
        strictver (bool): Strictly check the version requirements for
            this version of the dap.  If True, an exception is raised if
            the expected version are not the same as the environment.

    """
    def __init__(self, version=None, strictver=True):
        # Default to MPL-6
        self.mplver = 'MPL-6' if version is None else version

        mpls = self._available_mpls()
        mpli = numpy.where(mpls[:,0] == self.mplver)
        if len(mpli[0]) == 0:
            mpls = self._available_mpls(write=True)
            raise ValueError('{0} is not an available MPL!'.format(self.mplver))
           
        n = mpls.shape[0] 
        i = 0
        self.mplver = mpls[mpli,i].ravel()[0]
        i += 1
        self.accessver = mpls[mpli,i].ravel()[0]
        i += 1
        self.idlver = mpls[mpli,i].ravel()[0]
        i += 1
        self.corever = mpls[mpli,i].ravel()[0]
        i += 1
        self.drpver = mpls[mpli,i].ravel()[0]
        i += 1
        self.dapver = mpls[mpli,i].ravel()[0]
        i += 1
        self.pythonver = mpls[mpli,i].ravel()[0]
        i += 1
        self.numpyver = mpls[mpli,i].ravel()[0]
        i += 1
        self.scipyver = mpls[mpli,i].ravel()[0]
        i += 1
        self.matplotlibver = mpls[mpli,i].ravel()[0]
        i += 1
        self.astropyver = mpls[mpli,i].ravel()[0]
        i += 1
        self.pydlver = mpls[mpli,i].ravel()[0]
        i += 1

        # Treat versioning strictly
        self.strictver = strictver


    def _available_mpls(self, write=False):
        """
        Return a list of the available MPLs to analyze, and printing it if requested.

        .. todo:
            Have this read a config file!

        Args:
            write (bool): Write the list of available MPLs to the
            stdout.

        Returns:
            numpy.array: A string array with the MPL dependencies
        """
        #                             MPL SDSS_ACCESS  IDLUTILS      CORE       DRP       DAP
        mpl_def = numpy.array([ [ 'MPL-1',      None, 'v5_5_16', 'v1_0_0', 'v1_0_0',     None,
        #                            python     numpy     scipy, matplotlib  astropy  pydl
                                       None,     None,     None,       None,    None,     None ],
                                [ 'MPL-2',      None, 'v5_5_17', 'v1_0_0', 'v1_1_2',     None,
                                       None,     None,     None,       None,    None,     None ],
                                ['v1_2_0',      None, 'v5_5_19', 'v1_1_0', 'v1_2_0',     None,
                                       None,     None,     None,       None,    None,     None ],
                                [ 'MPL-3',      None, 'v5_5_22', 'v1_1_0', 'v1_3_3', 'v1_0_0',
                                       None,     None,     None,       None,    None,     None ],
                                [ 'MPL-4',   '0.0.0', 'v5_5_23', 'v1_2_0', 'v1_5_1',  '1.1.1',
                                       None,     None,     None,       None,    None,     None ],
                                [ 'MPL-5',   '0.2.0', 'v5_5_26', 'v1_3_1', 'v2_0_1',  '2.0.2',
                                    '3.5.1', '1.11.0', '0.17.1',       None,  '1.1.2', '0.5.0' ],
                                [ 'trunk',   '0.2.1', 'v5_5_29', 'v1_5_0', 'v2_2_0',  'trunk',
                                    '3.5.3', '1.12.0', '0.19.0',    '2.0.0',  '1.3.0', '0.5.3' ]
                              ])
        nmpl = mpl_def.shape[1]
        if write:
            print('Available MPLs:')
            print('{0:>7} {1:>7} {2:>8} {3:>6} {4:>6} {5:>6} {6:>6} {7:>6} {8:>6} {9:>10}'
                  ' {10:>7} {11:>6}'.format('MPL', 'SDSSACC', 'IDLUTILS', 'CORE', 'DRP', 'DAP',
                                            'PYTHON', 'NUMPY', 'SCIPY', 'MATPLOTLIB', 'ASTROPY',
                                            'PYDL'))
            for x in mpl_def[0:nmpl,:]:
                _x = x.astype(str)
                print('{0:>7} {1:>7} {2:>8} {3:>6} {4:>6} {5:>6} {6:>6} {7:>6} {8:>6} {9:>10}'
                      ' {10:>7} {11:>6}'.format(*_x))

        return mpl_def


    def _version_mismatch(self, package, required_ver, environment_ver, quiet=False):
        """
        A version mismatch was detected.  Handle based on
        :attr:`strictver`.
        """
        if self.strictver:
            raise EnvironmentError('Environment mismatch for {0}: expected {1}, found {2}'.format(
                                   package, required_ver, environment_ver))
        elif not quiet:
            warnings.warn('Environment mismatch for {0}: expected {1}, found {2}'.format(package,
                          required_ver, environment_ver))


    def verify_versions(self, quiet=True):
        """
        Verify that the defined MPL versions are consistent with the
        system where the DAP is being run.
        """
        # These will throw KeyErrors if the appropriate environmental variables do not exist
        # TODO: Use the python3 compatible versions of sdss4tools and sdss_access
        try:
            accessver_env = os.environ['SDSS_ACCESS_DIR'].split('/')[-1]
        except KeyError:
            accessver_env = None
        idlver_env = os.environ['IDLUTILS_DIR'].split('/')[-1]
        python_ver = '.'.join([ str(v) for v in sys.version_info[:3]])

        if not quiet:
            print('Current environment: ')
            print('    SDSS_ACCESS: {0}'.format(accessver_env))
            print('       IDLUTILS: {0}'.format(idlver_env))
            print('      MANGACORE: {0}'.format(os.environ['MANGACORE_VER']))
            print('       MANGADRP: {0}'.format(os.environ['MANGADRP_VER']))
            print('       MANGADAP: {0}'.format(os.environ['MANGADAP_VER']))
            print('         PYTHON: {0}'.format(python_ver))
            print('          NUMPY: {0}'.format(numpy.__version__))
            print('          SCIPY: {0}'.format(scipy.__version__))
            print('     MATPLOTLIB: {0}'.format(matplotlib.__version__))
            print('        ASTROPY: {0}'.format(astropy.__version__))
            print('           PYDL: {0}'.format(pydl.__version__))

        # Check versions
        if self.accessver != accessver_env:
            self._version_mismatch('SDSS_ACCESS', self.accessver, accessver_env, quiet=quiet)

        if self.idlver != idlver_env:
            self._version_mismatch('IDLUTILS', self.idlver, idlver_env, quiet=quiet)

        if self.corever != os.environ['MANGACORE_VER']:
            self._version_mismatch('MANGACORE', self.corever, os.environ['MANGACORE_VER'],
                                   quiet=quiet)

        if self.drpver != os.environ['MANGADRP_VER']:
            self._version_mismatch('MANGADRP', self.drpver, os.environ['MANGADRP_VER'], quiet=quiet)

        if self.pythonver is not None \
                and StrictVersion(python_ver) != StrictVersion(self.pythonver):
            self._version_mismatch('python', self.pythonver, python_ver, quiet=quiet)

        if self.numpyver is not None \
                and StrictVersion(numpy.__version__) != StrictVersion(self.numpyver):
            self._version_mismatch('numpy', self.numpyver, numpy.__version__, quiet=quiet)

        if self.scipyver is not None \
                and StrictVersion(scipy.__version__) != StrictVersion(self.scipyver):
            self._version_mismatch('scipy', self.scipyver, scipy.__version__, quiet=quiet)

        if self.astropyver is not None \
                and StrictVersion(astropy.__version__) != StrictVersion(self.astropyver):
            self._version_mismatch('astropy', self.astropyver, astropy.__version__, quiet=quiet)

        if self.pydlver is not None \
                and StrictVersion(pydl.__version__) != StrictVersion(self.pydlver):
            self._version_mismatch('pydl', self.pydlver, pydl.__version__, quiet=quiet)


#    def module_file(self):
#        """
#        Return the name of the module file specific to the MPL to analyze.
#
#        .. warning::
#            This function is incomplete in the sense that not all MPLs
#            listed in the above table have an associated module file.
#
#        Returns:
#            str: Name of the manga module to load.  E.g. 'manga/mpl3'
#            for MPL-3.
#        """
#        if self.mplver == 'MPL-1':
#            return 'manga/westfall3_mpl1'
#        elif self.mplver == 'MPL-2':
#            return 'manga/westfall3_mpl2'
#        elif self.mplver == 'MPL-3':
#            return 'manga/mpl3'
#        elif self.mplver == 'MPL-4':
#            return 'manga/mpl4'
#        elif self.mplver == 'MPL-5':
#            return 'manga/mpl5'
#        elif self.mplver == 'MPL-6':
#            return 'manga/mpl6'
#        else:
#            return 'manga/westfall3_mpl1'


    def show(self):
        """Print the available MPL versions to stdout."""

        print('{0}: SDSS_ACCESS:{1}; IDLUTILS:{2}; COREVER:{3}; DRPVER:{4}; DAPVER:{5}'.format(
              self.mplver, self.accessver, self.idlver, self.corever, self.drpver, self.dapver))


