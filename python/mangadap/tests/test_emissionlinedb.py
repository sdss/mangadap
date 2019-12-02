
import pytest
import os

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from mangadap.par.emissionlinedb import EmissionLineDB

def test_instantiation():

    emldb = EmissionLineDB.from_key('ELPMPL9')
    assert emldb.neml == 35

    _emldb = EmissionLineDB(emldb.file)
    assert _emldb.neml == 35

def test_channel_names():
    emldb = EmissionLineDB.from_key('ELPMPL9')
    names = emldb.channel_names()
    assert 'Ha-6564' in list(names.keys())
    assert names['OII-3727'] == 0

