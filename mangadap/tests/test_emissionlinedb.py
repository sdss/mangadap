from IPython import embed

import numpy

from mangadap.par import emissionlinedb
from mangadap.tests.util import data_test_file


def test_read():
    dbs = emissionlinedb.EmissionLineDB.available_databases()
    assert len(dbs) > 0, 'No emission-line databases available'
    for key in dbs.keys():
        emldb = emissionlinedb.EmissionLineDB.from_key(key)


def test_mpl11():
    emldb = emissionlinedb.EmissionLineDB.from_key('ELPMPL11')
    assert len(emldb) == 35, 'Incorrect number of emission lines'
    assert 'ArIII' in emldb['name'], 'Does not contain ArIII in list'


def test_datatable():

    emldb = emissionlinedb.EmissionLineDB(data_test_file('elp_ha_tied.par'))
    assert emldb['index'].shape == (5,), 'Incorrect number of lines'
    assert emldb['tie_par'].shape == (5,3), 'Incorrect number of tying parameters'
    assert numpy.array_equal(emldb['tie_index'][0],[3,3,3]) , \
        'All parameters of NII line should be tied to its doublet'

    indx_match = emldb.tie_index_match()
    assert numpy.array_equal(indx_match[0],[2,2,2]) , 'Index matching failed'

    tbl = emldb.to_datatable()
    assert numpy.array_equal(tbl['TIE_ID'], emldb['tie_index']), 'Bad conversion to datatable'

