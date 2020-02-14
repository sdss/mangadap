
import pytest

from mangadap.proc.spatiallybinnedspectra import available_spatial_binning_methods

def test_avail():
    methods_list = available_spatial_binning_methods()
    assert len(methods_list) > 0, 'No binning methods available'

