# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Classes to handle MaNGA DAP analysis plans.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
from pydl.pydlutils.yanny import yanny

from .parset import KeywordParSet, ParDatabase

class AnalysisPlan(KeywordParSet):
    """
    Generic class to handle MaNGA DAP analysis plans.

    The defined parameters are:

    .. include:: ../tables/analysisplan.rst

    """
    def __init__(self, drpqa_key=None, drpqa_clobber=None, bin_key=None, bin_clobber=None,
                 continuum_key=None, continuum_clobber=None, elmom_key=None, elmom_clobber=None,
                 elfit_key=None, elfit_clobber=None, spindex_key=None, spindex_clobber=None):

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
        defaults = [ None, False, None, False, None, False, None, False, None, False, None, False ]

        descr = ['Data reduction quality assessment method keyword',
                 'Overwrite any existing data-quality assessment reference files',
                 'Spatial binning method keyword',
                 'Overwrite any existing spatial binning reference files',
                 'Stellar-continuum fitting method keyword',
                 'Overwrite any existing stellar-continuum fitting reference files',
                 'Emission-line moments measurement method keyword',
                 'Overwrite any existing emission-line moments reference files',
                 'Emission-line modeling method keyword',
                 'Overwrite any existing emission-line modeling reference files',
                 'Spectral-index measurement method keyword',
                 'Overwrite any existing spectral-index reference files']

        super(AnalysisPlan, self).__init__(pars, values=values, defaults=defaults, dtypes=dtypes,
                                           descr=descr)
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

        super(AnalysisPlanSet, self).__init__(_planlist)

    @classmethod
    def from_par_file(cls, f):
        """
        Instantiate the plan set from a par file
        """

        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...

        if not os.path.isfile(f):
            raise FileNotFoundError('No file {0}.'.format(f))
    
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

    @classmethod
    def default(cls):
        """
        Return the default analysis plan set.
        """
        return cls([AnalysisPlan(drpqa_key='SNRG', bin_key='HYB10', continuum_key='MILESHCMPL11',
                                 elmom_key='EMOMMPL11', elfit_key='EFITMPL11HCDB',
                                 spindex_key='INDXEN')])

