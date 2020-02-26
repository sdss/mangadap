from IPython import embed
import pytest

from mangadap.proc.stellarcontinuummodel import available_stellar_continuum_modeling_methods
from mangadap.proc.templatelibrary import TemplateLibrary

def test_avail():
    methods_list = available_stellar_continuum_modeling_methods()
    assert len(methods_list) > 0, 'No stellar-continuum modeling methods available'

    for method in methods_list:
        if 'template_library' in method['fitpar'].keys():
            tpl = TemplateLibrary(method['fitpar']['template_library_key'], match_resolution=False,
                                  velscale_ratio=1, spectral_step=1e-4, log=True, hardcopy=False)

