# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a few base classes used during spectral fitting procedures.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spectralfitting.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        import warnings
        if sys.version > '3':
            long = int

        import numpy

*Class usage examples*:
        Add examples

*Revision history*:
    | **14 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **19 Apr 2016**: (KBW) First version
    | **26 Apr 2016**: (KBW) Moved PPXFFit to a separate file (ppxffit.py)

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import numpy
from ..util.bitmask import BitMask

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion



# BASE CLASS -----------------------------------------------------------

class SpectralFitting():
    """
    Base class for spectral fitting.
    """
    def __init__(self, fit_type, bitmask, par=None):
        self.fit_type = fit_type
        if not isinstance(bitmask, BitMask):
            raise TypeError('Input bit mask must have type BitMask.')
        self.bitmask = bitmask
        self.par = par

# ----------------------------------------------------------------------
class StellarKinematicsFit(SpectralFitting):
    """
    Base class for fitting stellar kinematics.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'stellar_kinematics', bitmask, par=par)
        self.fit_method = fit_method

    @staticmethod
    def _per_stellar_kinematics_dtype(ntpl, nadd, nmult, nkin, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('MASK',mask_dtype),
                 ('BEGPIX', numpy.int),
                 ('ENDPIX', numpy.int),
                 ('NPIXTOT',numpy.int),
                 ('NPIXFIT',numpy.int),
                 ('TPLWGT',numpy.float,ntpl),
                 ('ADDCOEF',numpy.float,nadd) if nadd > 1 else ('ADDCOEF',numpy.float),
                 ('MULTCOEF',numpy.float,nmult) if nmult > 1 else ('MULTCOEF',numpy.float),
                 ('KININP',numpy.float,2),
                 ('KIN',numpy.float,nkin),
                 ('KINERR',numpy.float,nkin),
                 ('CHI2',numpy.float),
                 ('RCHI2',numpy.float),
                 ('ROBUST_RCHI2',numpy.float),
                 ('RMS',numpy.float),
                 ('ABSRESID',numpy.float,5),
                 ('FRMS',numpy.float),
                 ('FABSRESID',numpy.float,5),
                 ('SIGMACORR',numpy.float)
               ]

# ----------------------------------------------------------------------
class CompositionFit(SpectralFitting):
    """
    Base class for fitting the spectral composition.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'composition', bitmask, par=par)
        self.fit_method = fit_method


# ----------------------------------------------------------------------
class EmissionLineFit(SpectralFitting):
    """
    Base class for fitting emission lines.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'emission_line', bitmask, par=par)
        self.fit_method = fit_method


    @staticmethod
    def _per_emission_line_dtype(neml, nkin, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('WIN_INDEX',numpy.int,neml),
                 ('MASK', mask_dtype,neml),
                 ('FLUX',numpy.float,neml),
                 ('FLUXERR',numpy.float,neml),
                 ('KIN',numpy.float,(neml,nkin)),
                 ('KINERR',numpy.float,(neml,nkin)),
                 ('SINST',numpy.float,neml)
               ]

                 
#                 ('EWOF',numpy.float,neml),
#                 ('EWOFERR',numpy.float,neml)
#                 ('EW',numpy.float,neml),
#                 ('EWERR',numpy.float,neml)
#                 ('COMBKIN',numpy.float,nkin),
#                 ('COMBKINERR',numpy.float,nkin),
#                 ('COMBKINSTD',numpy.float,nkin),
#                 ('COMBKINN',numpy.int),

    @staticmethod
    def _per_fitting_window_dtype(nwin, max_npar, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """

        return [ ('BIN_INDEX',numpy.int),
                 ('MASK', mask_dtype, nwin),
                 ('NPIXFIT',numpy.int,nwin),
                 ('PAR',numpy.float,(nwin,max_npar)),
                 ('ERR',numpy.float,(nwin,max_npar)),
                 ('LOBND',numpy.float,(nwin,max_npar)),
                 ('UPBND',numpy.float,(nwin,max_npar)),
                 ('FIXED',numpy.bool,(nwin,max_npar)),
#                 ('TIED',numpy.int,(nwin,max_npar)),
                 ('IGNORE',numpy.bool,(nwin,max_npar)),
                 ('CHI2',numpy.float,nwin),
                 ('RCHI2',numpy.float,nwin),
                 ('RMS',numpy.float,nwin),
                 ('RESID',numpy.float,(nwin,7)),
                 ('FRAC_RMS',numpy.float,nwin),
                 ('FRAC_RESID',numpy.float,(nwin,7))
               ]


