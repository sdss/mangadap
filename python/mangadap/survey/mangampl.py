"""
Defines a class to keep track of MaNGA MPL dependencies and versioning.

.. note::
    This class was created for the survey-level execution of the DAP.
    The dependencies will necessarily grow, and will be handled using
    module files.  So, even as it is, this class has little utility.

*License*:
    Copyright (c) 2015, Kyle B. Westfall, Brett H. Andrews
    Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/mangampl.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import numpy

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **2015**: Original Implementation by K. Westfall (KBW)
    | **Apr 2015**: (KBW) Added MPL-3
    | **20 May 2015**: (KBW) Documentation and Sphinx tests. Change
        default to MPL-3
    | **27 Aug 2015**: (KBW) Added (temporary) MPL-4 dependencies
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy

# TODO: Add DAP version

class mangampl:
    """
    Create an object that keeps track of IDLUTILS and DRP dependencies
    for a given MPL and other intermediate analyses.  Dependencies are
    selected by the "MPL" version as provided in the following table:

    +--------+------------+--------+
    |  MPL   |  IDLUTILS  |   DRP  |
    +========+============+========+
    | MPL-1  | v5_5_16    | v1_0_0 |
    +--------+------------+--------+
    | MPL-2  | v5_5_17    | v1_1_2 |
    +--------+------------+--------+
    | v1_2_0 | v5_5_19    | v1_2_0 |
    +--------+------------+--------+
    | MPL-3  | v5_5_22    | v1_3_3 |
    +--------+------------+--------+
    | MPL-4  | v5_5_22    | v1_4_0 |
    +--------+------------+--------+

    .. warning::
        Definition of MPL-4 is only temporary!

    Args:
        version (str): Version as selected via matching to the "MPL"
            string in the table listed above. Default: MPL-3.

    Raises:
        Exception: Raised if the provided *version* is not one of the
            MPL values in the above table.

    Attributes:
        mplver (str): Version as selected via matching to the "MPL"
            string in the table listed above. Default: MPL-3.
        idlver (str): IDLUTILS version.
        drpver (str): MANGADRP version.

    """
    def __init__(self, version=None):

        # Default to MPL-4
        self.mplver = 'MPL-4' if version is None else version

        mpls = self._available_mpls()
        mpli = numpy.where(mpls[:,0] == self.mplver)
        if len(mpli[0]) == 0:
            mpls = self._available_mpls(write=True)
            raise Exception('{0} is not an available MPL!'.format(self.mplver))
            
        self.mplver = str(mpls[mpli].reshape(3)[0])
        self.idlver = str(mpls[mpli].reshape(3)[1])
        self.drpver = str(mpls[mpli].reshape(3)[2])


    def _available_mpls(self, write=False):
        """
        Return a list of the available MPLs to analyze, and printing it if requested.

        Args:
            write (bool): Write the list of available MPLs to the
            stdout.

        Returns:
            numpy.array : A string array with the MPL dependencies

        """
        nmpl = 5
        #                        MPL      IDLUTILS   DRP
        mpl_def = numpy.array([
                                ['MPL-1', 'v5_5_16', 'v1_0_0'],
                                ['MPL-2', 'v5_5_17', 'v1_1_2'],
                                ['v1_2_0', 'v5_5_19', 'v1_2_0'],
                                ['MPL-3', 'v5_5_22', 'v1_3_3'],
                                ['MPL-4', 'v5_5_22', 'v1_4_0']
                              ])

        if write:
            print('Available MPLs:')
            for x in mpl_def[0:nmpl,:]:
                print('{0}: IDLUTILS:{1}; DRPVER:{2}'.format(x[0], x[1], x[2]))

        return mpl_def


    def module_file(self):
        """
        Return the name of the module file specific to the MPL to analyze.

        .. warning::
            This function is incomplete in the sense that not all MPLs
            listed in the above table have an associated module file.

        Returns:
            str : Name of the manga module to load.  E.g. 'manga/mpl3'
                for MPL-3.

        """
        if self.mplver == 'MPL-1':
            return 'manga/westfall3_mpl1'
        elif self.mplver == 'MPL-2':
            return 'manga/westfall3_mpl2'
        elif self.mplver == 'MPL-3':
            return 'manga/mpl3'
        elif self.mplver == 'MPL-4':
            return 'manga/mpl4'
        else:
            return 'manga/westfall3_mpl1'


    def show(self):
        """Print the available MPL versions to stdout."""

        print('{0}: IDLUTILS:{1}; DRPVER:{2}'.format(self.mplver, self.idlver, self.drpver))



