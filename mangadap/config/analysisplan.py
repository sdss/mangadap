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

from ..par.parset import KeywordParSet, ParDatabase
from . import defaults
from ..proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from ..proc.stellarcontinuummodel import StellarContinuumModel
from ..proc.emissionlinemodel import EmissionLineModel

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

        key = AnalysisPlan.unique_plan_key(_bin_key, _continuum_key, _elfit_key)

        pars =   ['key', 'drpqa_key', 'drpqa_clobber', 'bin_key', 'bin_clobber', 'continuum_key',
                  'continuum_clobber', 'elmom_key', 'elmom_clobber', 'elfit_key', 'elfit_clobber',
                  'spindex_key', 'spindex_clobber']
        values = [key, _drpqa_key, drpqa_clobber, _bin_key, bin_clobber, _continuum_key,
                  continuum_clobber, _elmom_key, elmom_clobber, _elfit_key, elfit_clobber,
                  _spindex_key, spindex_clobber]
        dtypes = [str, str, bool, str, bool, str, bool, str, bool, str, bool, str, bool]
        defaults = [None, None, False, None, False, None, False, None, False, None, False, None,
                    False]
        descr = ['Unique key used to identify this plan',
                 'Data reduction quality assessment method keyword',
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

    @staticmethod
    def unique_plan_key(bin_key, continuum_key, elfit_key):
        bin_method = SpatiallyBinnedSpectra.define_method(bin_key)
        sc_method = StellarContinuumModel.define_method(continuum_key)
        el_method = EmissionLineModel.define_method(elfit_key)
        eltpl = sc_method["fitpar"]["template_library_key"] \
                    if el_method['continuum_tpl_key'] is None else el_method['continuum_tpl_key']
        return f'{bin_method["key"]}-{sc_method["fitpar"]["template_library_key"]}-{eltpl}'

        
class AnalysisPlanSet(ParDatabase):
    """
    Container class for a set of analysis plans.
    """
    def __init__(self, planlist, cube=None, analysis_path=None):
        _planlist = planlist if isinstance(planlist, list) else [planlist]
        self.nplans = len(planlist)
        for i in range(self.nplans):
            if not isinstance(_planlist[i], AnalysisPlan):
                raise TypeError('Input must be a single or list of AnalysisPlan objects.')
        super().__init__(_planlist)
        self.root_path = Path('.' if analysis_path is None else analysis_path).resolve()
        self.cube = cube

    @classmethod
    def from_par_file(cls, f, **kwargs):
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
        planlist = [AnalysisPlan(par['DAPPLAN']['drpqa_key'][i],
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
                                 bool(par['DAPPLAN']['spindex_clobber'][i]))
                        for i in range(nplan)]
        return cls(planlist, **kwargs)

    @classmethod
    def default(cls, **kwargs):
        """
        Return the default analysis plan set.
        """
        return cls([AnalysisPlan(drpqa_key='SNRG', bin_key='HYB10', continuum_key='MILESHCMPL11',
                                 elmom_key='EMOMMPL11', elfit_key='EFITMPL11HCDB',
                                 spindex_key='INDXEN')], **kwargs)

    def common_path(self):
        """
        Return the path for data common to all plans.

        Args:
            cube (:class:`~mangadap.datacube.datacube.DataCube`, optional):
                Cube being analyzed.  Passed for cube-specific path
                specification.  Not used by this base class.

        Returns:
            `Path`_: Path object for the "common" output
        """
        return self.root_path / 'common'

    def method_path(self, plan_index=0, qa=False, ref=False):
        """
        Return the path for method-specific output.

        Args:
            cube (:class:`~mangadap.datacube.datacube.DataCube`, optional):
                Cube being analyzed.  Passed for cube-specific path
                specification.  Not used by this base class.
            plan_index (:obj:`int`, optional):
                The index of the plan.  This is used to select the 'key' of the
                analysis plan being used, which is used as the subdirectory for
                the output.
            qa (:obj:`bool`, optional):
                Flag that the output is specifically quality assessment plots.
            ref (:obj:`bool`, optional):
                Flag that the output is specifically reference data.

        Returns:
            `Path`_: Path object for the method-specific output

        Raises:
            ValueError:
                Raised if the plan index is invalid or if both qa and ref are true.

        """
        if plan_index < 0 or plan_index >= self.nplan:
            raise ValueError(f'Invalid index ({plan_index}); must be 0 < index < {self.nplan}.')
        if qa and ref:
            raise ValueError('Cannot provide path for both qa and ref directory.  Pick one.')
        root = self.root_path / self['key'][plan_index]
        if not qa and not ref:
            return root
        if qa:
            return root / 'qa'
        return root / 'ref'

    def dap_file_root(self, cube, mode=None, plan_index=None):
        """
        Return the general root name for an output file.

        The default returned by this base class is::

        {cube.output_root}-{mode}-{self['key'][plan_index]}
        
        where ``mode`` and ``self['key'][plan_index]`` are only included if
        ``mode`` or ``plan_index`` are provided.

        Args:
            cube (:class:`~mangadap.datacube.datacube.DataCube`):
                Cube being analyzed.  This is used to provide the primary root
                of the output name.
            mode (:obj:`str`, optional):
                An optional "mode" included in the file name.  If None, not
                include in the output string.
            plan_index (:obj:`int`, optional):
                The index of the plan.  This is used to select the 'key' of the
                analysis plan being used, which is included in the file name.
                If None, this is not included in the file name.
    
        Returns:
            :obj:`str`: General root for output DAP file names.
        """
        root = cube.output_root
        if mode is not None:
            root = f'{root}-{mode}'
        if plan_index is not None:
            if plan_index < 0 or plan_index >= self.nplan:
                raise ValueError(f'Invalid index ({plan_index}); must be 0 < index < {self.nplan}.')
            root = f'{root}-{self["key"][plan_index]}'
        return root


