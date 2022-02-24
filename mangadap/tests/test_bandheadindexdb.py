
from IPython import embed

from mangadap.par.bandheadindexdb import BandheadIndexDB


def test_read():
    dbs = BandheadIndexDB.available_databases()
    for key in dbs.keys():
        bhddb = BandheadIndexDB.from_key(key)


def test_basic():
    bhddb = BandheadIndexDB.from_key('BHBASIC')
    assert len(bhddb) == 3, 'Incorrect number of bandhead indices'
    assert 'Dn4000' in bhddb['name'], 'Does not contain Dn4000 in list'


