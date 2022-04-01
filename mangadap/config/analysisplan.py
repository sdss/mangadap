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

from pathlib import Path
from copy import deepcopy

from IPython import embed

import tomli

from ..par.util import recursive_dict_str_to_None
from . import defaults

from mangadap.proc.reductionassessments import ReductionAssessmentDef
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectraDef
from mangadap.proc.stellarcontinuummodel import StellarContinuumModelDef
from mangadap.proc.emissionlinemoments import EmissionLineMomentsDef
from mangadap.proc.emissionlinemodel import EmissionLineModelDef
from mangadap.proc.spectralindices import SpectralIndicesDef


class AnalysisPlan:
    """
    Container class for a set of analysis plans.
    """
    def __init__(self, plan, cube=None, analysis_path=None):
        if not isinstance(plan, dict):
            raise TypeError('Plan must be provided as a dictionary.')
        # Copy the provided plan dictionary so that it can be changed/filled
        # throughout the execution of the code.
        self.plan = deepcopy(plan)
        self.plan_keys = list(self.plan.keys())
        self.nplans = len(self.plan_keys)
        # Make sure that the plans have a key name.  If not, just use the
        # dictionary key of the dictionary for each plan.
        for key in self.plan_keys:
            if 'key' not in self.plan[key]:
                self.plan['key'] = key
        self.analysis_path = Path('.' if analysis_path is None else analysis_path).resolve()
        self.cube = cube
        self._validate()
        self.parse()

    def _validate(self):
        """
        Validate the provided plans.
        """
        # For now, just make sure that all the plans have the base top-level keys
        required_keys = ['rdxqa', 'binning', 'continuum', 'eline_moments', 'eline_fits', 'indices']
        for key in self.plan_keys:
            for method_key in required_keys:
                if method_key not in self[key].keys():
                    self.plan[key][method_key] = {}

        # Recursively convert None strings into None types
        self.plan = recursive_dict_str_to_None(self.plan)

        # TODO: Cross validate parameter sets between plans!!

    def __getitem__(self, key):
        """
        Return the value of the designated key.
        """
        return self.plan[key]

    def keys(self):
        return self.plan.keys()

    @classmethod
    def from_toml(cls, ifile, **kwargs):
        """
        Instantiate the plan from a TOML file.
        """
        _ifile = Path(ifile).resolve()
        if not _ifile.exists():
            raise FileNotFoundError(f'{_ifile} does not exist!')

        with open(_ifile, 'rb') as f:
            plan = tomli.load(f)

        return cls(plan, **kwargs)

    @classmethod
    def default(cls, **kwargs):
        """
        Return the default analysis plan set.
        """
        f = defaults.dap_config_root() / 'default_plan.toml'
        return cls.from_toml(f, **kwargs)

    def parse(self):
        self.rdxqa = {key: None if self[key]['rdxqa'] is None else
                        ReductionAssessmentDef.from_dict(self[key]['rdxqa'])
                        for key in self.plan.keys()}
        self.binning = {key: None if self[key]['binning'] is None else
                            SpatiallyBinnedSpectraDef.from_dict(self[key]['binning'])
                            for key in self.plan.keys()}
        self.continuum = {key: None if self[key]['continuum'] is None else
                            StellarContinuumModelDef.from_dict(self[key]['continuum'])
                            for key in self.plan.keys()}
        self.elmom = {key: None if self[key]['eline_moments'] is None else
                        EmissionLineMomentsDef.from_dict(self[key]['eline_moments']) 
                        for key in self.plan.keys()}
        self.elfit = {key: None if self[key]['eline_fits'] is None else
                        EmissionLineModelDef.from_dict(self[key]['eline_fits'])
                        for key in self.plan.keys()}
        self.sindx = {key: None if self[key]['indices'] is None else
                        SpectralIndicesDef.from_dict(self[key]['indices'])
                        for key in self.plan.keys()}

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
        return self.analysis_path / 'common'

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
        if plan_index < 0 or plan_index >= self.nplans:
            raise ValueError(f'Invalid index ({plan_index}); 0 <= index < {self.nplans}.')
        if qa and ref:
            raise ValueError('Cannot provide path for both qa and ref directory.  Pick one.')
        root = self.analysis_path / self.plan[self.plan_keys[plan_index]]['key']
        if not qa and not ref:
            return root
        if qa:
            return root / 'qa'
        return root / 'ref'

    # TODO: Restrict mode?
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
            if plan_index < 0 or plan_index >= self.nplans:
                raise ValueError(f'Invalid index ({plan_index}): 0 <= index < {self.nplans}.')
            plan_name = list(self.keys())[plan_index]
            root = f'{root}-{self[plan_name]["key"]}'
        return root




