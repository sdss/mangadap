
import pytest
import os
from IPython import embed

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from mangadap.proc.emissionlinemodel import available_emission_line_modeling_methods

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.templatelibrary import TemplateLibrary

def test_avail():
    # Get the available methods
    method_list = available_emission_line_modeling_methods()
    assert len(method_list) > 0, 'No emission-line-modeling methods available'

    # For the available methods, make sure that the ancillary databases
    # and templates can be loaded.
    for method in method_list:
        if method['artifacts'] is not None:
            artdb = ArtifactDB.from_key(method['artifacts'])
        if method['ism_mask'] is not None:
            emldb = EmissionLineDB.from_key(method['ism_mask'])
        if method['emission_lines'] is not None:
            emldb = EmissionLineDB.from_key(method['emission_lines'])
        if method['continuum_tpl_key'] is not None:
            tpl = TemplateLibrary(method['continuum_tpl_key'], match_resolution=False,
                                  velscale_ratio=1, spectral_step=1e-4, log=True, hardcopy=False)
        if method['fitpar'] is None:
            continue

        if 'continuum_templates' in method['fitpar'].keys() \
                and method['fitpar']['continuum_templates'] is not None:
            tpl = TemplateLibrary(method['fitpar']['continuum_templates'],
                                  match_resolution=False, velscale_ratio=1,
                                  spectral_step=1e-4, log=True, hardcopy=False)
