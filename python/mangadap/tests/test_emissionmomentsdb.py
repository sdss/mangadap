
import pytest
import os

from IPython import embed

from mangadap.par.emissionmomentsdb import EmissionMomentsDB

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_read():
    dbs = EmissionMomentsDB.available_databases()
    assert len(dbs) > 0, 'No emission-line-moment databases found'
    for key in dbs.keys():
        print(key)
        emomdb = EmissionMomentsDB.from_key(key)

def test_mpl9():
    emomdb = EmissionMomentsDB.from_key('ELBMPL9')
    assert len(emomdb) == 35, 'Incorrect number of emission-line bands'
    assert 'ArIII' in emomdb['name'], 'Does not contain ArIII in list'

