# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class to keep track of MaNGA MPL dependencies and versioning.

.. note::
    This class was created for the survey-level execution of the DAP.
    The dependencies will necessarily grow, and will be handled using
    module files.  So, even as it is, this class has little utility.

.. todo::
    Deprecate this. I'm not sure it's actually used (or useful)
    anymore.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import sys
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
    | MPL-6  |        None |    v5_5_30 | v1_6_0 | v2_3_1 |  2.1.3 |
    +--------+-------------+------------+--------+--------+--------+
    | MPL-7  |        None |    v5_5_32 | v1_6_2 | v2_4_3 |  2.2.1 |
    +--------+-------------+------------+--------+--------+--------+

    **Python versions**
    +--------+--------+--------+--------+------------+---------+-------+----------------+
    | MPL    | python | numpy  |  scipy | matplotlib | astropy |  pydl | Fwd compatible |
    +========+========+========+========+============+=========+=======+================+
    | MPL-5  |  3.5.1 | 1.11.0 | 0.17.1 |       None |   1.1.2 | 0.5.0 |            yes |
    +--------+--------+--------+--------+------------+---------+-------+----------------+
    | MPL-6  |  3.6.3 | 1.13.3 |  1.0.0 |       None |   2.0.2 | 0.6.0 |            N/A |
    +--------+--------+--------+--------+------------+---------+-------+----------------+
    | MPL-7  |  3.6.3 | 1.13.3 |  1.0.0 |       None |   2.0.2 | 0.6.0 |            N/A |
    +--------+--------+--------+--------+------------+---------+-------+----------------+

    .. note::
        - "Fwd compatible" means that you can run this version of the
          DAP with the most recent version without needed to revert the
          python packages to their previous versions.
        - Only the survey-level routines are dependent on MANGACORE; the
          core DAP processing software only depends on the listed python
          packages.

    .. todo::
        - Change python dependencies to "minimum" versions?

    .. warning::
        - SDSS_ACCESS dependencies are only because of the plotting
          package that has since moved into Marvin.
        - For MPL-1, the MANGACORE version was set to trunk.  Here it
          defaults to v1_0_0.
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
        # Default to current environment
        self.mplver = 'env' if version is None else version

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
        #                            python     numpy     scipy  matplotlib  astropy      pydl
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
                                [ 'MPL-6',      None, 'v5_5_30', 'v1_6_0', 'v2_3_1',  '2.1.3',
                                    '3.6.3', '1.13.3',    '1.0',       None,  '2.0.2', '0.6.0' ],
                                [ 'MPL-7',      None, 'v5_5_32', 'v1_6_2', 'v2_4_3',  '2.2.1',
                                    '3.6.3', '1.13.3',    '1.0',       None,  '2.0.2', '0.6.0' ],
                                ['env'] + self.get_environment_versions()
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


    def get_environment_versions(self):
        # These will throw KeyErrors if the appropriate environmental variables do not exist
        try:
            accessver_env = os.environ['SDSS_ACCESS_DIR'].split('/')[-1]
        except KeyError:
            accessver_env = None
        try:
            idlver_env = os.environ['IDLUTILS_DIR'].split('/')[-1]
        except KeyError:
            idlver_env = None
        python_ver = '.'.join([ str(v) for v in sys.version_info[:3]])

        return [ accessver_env, idlver_env, os.getenv('MANGACORE_VER'),
                 os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'),
                python_ver, numpy.__version__, scipy.__version__, matplotlib.__version__,
                astropy.__version__, pydl.__version__]


    def verify_versions(self, quiet=True):
        """
        Verify that the defined MPL versions are consistent with the
        system where the DAP is being run.
        """
        env = self.get_environment_versions()

        if not quiet:
            print('Current environment: ')
#            print('    SDSS_ACCESS: {0}'.format(env[0]))
            print('       IDLUTILS: {0}'.format(env[1]))
            print('      MANGACORE: {0}'.format(env[2]))
            print('       MANGADRP: {0}'.format(env[3]))
            print('       MANGADAP: {0}'.format(env[4]))
            print('         PYTHON: {0}'.format(env[5]))
            print('          NUMPY: {0}'.format(env[6]))
            print('          SCIPY: {0}'.format(env[7]))
            print('     MATPLOTLIB: {0}'.format(env[8]))
            print('        ASTROPY: {0}'.format(env[9]))
            print('           PYDL: {0}'.format(env[10]))

        # Check versions
#        if self.accessver != env[0]:
#            self._version_mismatch('SDSS_ACCESS', self.accessver, env[0], quiet=quiet)

        if self.idlver != env[1]:
            self._version_mismatch('IDLUTILS', self.idlver, env[1], quiet=quiet)

        if self.corever != env[2]:
            self._version_mismatch('MANGACORE', self.corever, env[2], quiet=quiet)

        if self.drpver != env[3]:
            self._version_mismatch('MANGADRP', self.drpver, env[3], quiet=quiet)

#        if self.dapver != env[4]:
#            self._version_mismatch('MANGADAP', self.dapver, env[4], quiet=quiet)

        if self.pythonver is not None and StrictVersion(env[5]) != StrictVersion(self.pythonver):
            self._version_mismatch('python', self.pythonver, env[5], quiet=quiet)

        if self.numpyver is not None and StrictVersion(env[6]) != StrictVersion(self.numpyver):
            self._version_mismatch('numpy', self.numpyver, env[6], quiet=quiet)

        if self.scipyver is not None and StrictVersion(env[7]) != StrictVersion(self.scipyver):
            self._version_mismatch('scipy', self.scipyver, env[7], quiet=quiet)

#        if self.matplotlibver is not None \
#                and StrictVersion(env[8]) != StrictVersion(self.matplotlibver):
#            self._version_mismatch('matplotlib', self.matplotlibver, env[8], quiet=quiet)

        if self.astropyver is not None and StrictVersion(env[9]) != StrictVersion(self.astropyver):
            self._version_mismatch('astropy', self.astropyver, env[9], quiet=quiet)

        if self.pydlver is not None and StrictVersion(env[10]) != StrictVersion(self.pydlver):
            self._version_mismatch('pydl', self.pydlver, env[10], quiet=quiet)


    def show(self):
        """Print the available MPL versions to stdout."""

        print('{0}: IDLUTILS:{1}; COREVER:{2}; DRPVER:{3}; DAPVER:{4}'.format(
              self.mplver, self.idlver, self.corever, self.drpver, self.dapver))


