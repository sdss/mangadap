# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Classes to handle MaNGA DAP analysis plans.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/analysisplan.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import numpy
        from pydl.pydlutils.yanny import yanny

        from .parset import ParSet, ParDatabase

*Class usage examples*:
        Add usage comments!

*Revision history*:
    | **11 May 2016**: Original implementation by K. Westfall (KBW)
    | **03 Jun 2016**: (KBW) Any key with value 'None' is converted to
        NoneType
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from pydl.pydlutils.yanny import yanny

from .parset import ParSet, ParDatabase

__author__ = 'Kyle B. Westfall'

class AnalysisPlan(ParSet):
    """
    Generic class to handle MaNGA DAP analysis plans.
    """

    def __init__(self, drpqa_key=None, drpqa_clobber=False, bin_key=None, bin_clobber=False,
                 continuum_key=None, continuum_clobber=False, elmom_key=None, elmom_clobber=False,
                 elfit_key=None, elfit_clobber=False, spindex_key=None, spindex_clobber=False):

        # TODO: Include covariance keys that is applied to each/some analyses
        _drpqa_key = None if drpqa_key == 'None' or drpqa_key is None else drpqa_key
        _bin_key = None if bin_key == 'None' or bin_key is None else bin_key
        _continuum_key = None if continuum_key == 'None' or continuum_key is None else continuum_key
        _elmom_key = None if elmom_key == 'None' or elmom_key is None else elmom_key
        _elfit_key = None if elfit_key == 'None' or elfit_key is None else elfit_key
        _spindex_key = None if spindex_key == 'None' or spindex_key is None else spindex_key

        pars =   [ 'drpqa_key', 'drpqa_clobber', 'bin_key', 'bin_clobber', 'continuum_key',
                   'continuum_clobber', 'elmom_key', 'elmom_clobber', 'elfit_key', 'elfit_clobber',
                   'spindex_key', 'spindex_clobber' ]
        values = [ _drpqa_key, drpqa_clobber, _bin_key, bin_clobber, _continuum_key,
                   continuum_clobber, _elmom_key, elmom_clobber, _elfit_key, elfit_clobber,
                   _spindex_key, spindex_clobber ]
        dtypes = [ str, bool, str, bool, str, bool, str, bool, str, bool, str, bool ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)
        self._check()


    def _check(self):
        """
        Check that the plan makes sense:
            - To bin, must first assess the DRP data
            - To do anythin else, must bin.
            - Emission-line moments and profile fits *should* subtract
              the best-fitting stellar continuum to account for the
              stellar-absorption below the lines, but it's not
              explicitly necessary for the code to run.
            - The spectral indices *should* subtract the best-fitting
              emission-line model and use the best-fitting stellar
              continuum to determine the velocity dispersion
              corrections, but neither is explicitly necessary for the
              code to run.
        """
        if self.data['bin_key'] is not None and self.data['drpqa_key'] is None:
            raise ValueError('To bin, must provide key for reduction assessments.')
        if self.data['continuum_key'] is not None and self.data['bin_key'] is None:
            raise ValueError('To fit the stellar continuum, must provide a binning key.')
        if self.data['elmom_key'] is not None and self.data['bin_key'] is None:
            raise ValueError('To measure the emission-line moments, must provide a binning key.')
        if self.data['elfit_key'] is not None and self.data['bin_key'] is None:
            raise ValueError('To measure the emission-line moments, must provide a binning key.')
        if self.data['spindex_key'] is not None and self.data['bin_key'] is None:
            raise ValueError('To measure the spectral indices, must provide a binning key.')

        
class AnalysisPlanSet(ParDatabase):
    """
    Container class for a set of analysis plans.
    """
    def __init__(self, planlist):
        _planlist = planlist if isinstance(planlist, list) else [planlist]
        self.nplans = len(planlist)
        for i in range(self.nplans):
            if not isinstance(_planlist[i], AnalysisPlan):
                raise TypeError('Input must be a single or list of AnalysisPlan objects.')

        ParDatabase.__init__(self, _planlist)

    @classmethod
    def from_par_file(cls, f):
        """
        Instantiate the plan set from a par file
        """

        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...
    
        par = yanny(filename=f, raw=True)
        if len(par['DAPPLAN']['drpqa_key']) == 0:
            raise ValueError('Could not find DAPPLAN entries in {0}!'.format(f))

        # Setup the array of emission line database parameters
        nplan = len(par['DAPPLAN']['drpqa_key'])
        planlist = []
        for i in range(nplan):
            planlist += [ AnalysisPlan(par['DAPPLAN']['drpqa_key'][i],
                                       bool(par['DAPPLAN']['drpqa_clobber'][i]),
                                       par['DAPPLAN']['bin_key'][i],
                                       bool(par['DAPPLAN']['bin_clobber'][i]),
                                       par['DAPPLAN']['continuum_key'][i],
                                       bool(par['DAPPLAN']['continuum_clobber'][i]),
                                       par['DAPPLAN']['elmom_key'][i],
                                       bool(par['DAPPLAN']['elmom_clobber'][i]),
                                       par['DAPPLAN']['elfit_key'][i],
                                       bool(par['DAPPLAN']['elfit_clobber'][i]),
                                       par['DAPPLAN']['spindex_key'][i],
                                       bool(par['DAPPLAN']['spindex_clobber'][i])) ]
        return cls(planlist)


#    SignalToNoisePlan
#        - wavelength range(s)
#        - covariance
#        - clobber
#    BinningPlan
#        - type: All
#        - type: None
#        - type: Voronoi
#            - S/N target
#            - S/N calculation
#        - type: Radial
#            - geometry: cx, cy, pa, boa
#            - rscale
#            - sampling: rs (in scale units), re (in scale units), nr, rlog
#        - velocity register
#        - weighting
#        - S/N Threshold
#        - rejection
#        - covariance
#        - clobber
#    StellarContinuumPlan: can be NULL
#        - (binned) spectra
#        - S/N Threshold
#        - template library
#        - pPXFPlan
#            - moments
#            - degree
#            - mdegree
#            - reddening
#            - oversample
#            ...
#        - clobber
#    EmissionLinePlan: can be NULL
#        - (binned) spectra
#        - S/N Threshold
#        - emission-line moment database
#        - emission-line database
#        - integrated flux limit
#        - allow negative amplitude?
#        - clobber
#    SpectralIndexPlan: can be NULL
#        - (binned) spectra
#        - S/N Threshold
#        - absorption-line index database
#        - bandhead index database
#        - pseudo-continuum limit?
#        - clobber
#
#    - index
#    - LIN vs LOG
#    - CUBE vs. RSS
#    - prior
#    - execute
#
#
#    open empty, open from file
#    add new plan, and iteratively do so
#    delete (all) plan(s)
#    print plans to screen or log file
#    write plans
#
