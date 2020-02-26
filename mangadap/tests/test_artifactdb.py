
import pytest
import os

from IPython import embed

from mangadap.par.artifactdb import ArtifactDB

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_read():
    dbs = ArtifactDB.available_databases()
    assert len(dbs) > 0, 'No artifact databases available'
    for key in dbs.keys():
        artdb = ArtifactDB.from_key(key)

def test_badsky():
    artdb = ArtifactDB.from_key('BADSKY')
    assert len(artdb) == 1, 'Incorrect number of bad sky bands'
    assert 'SKY' in artdb['name'], 'Does not contain SKY in list'

