
import pytest

from IPython import embed

import numpy
from astropy.io import fits

from mangadap.util.datatable import DataTable
from mangadap.proc.spectralfitting import StellarKinematicsFitDataTable
from mangadap.proc.spectralfitting import EmissionLineFitDataTable

def test_empty():
    keys = ['TEST', 'THIS']
    types = [ int, float]
    descr = ['This is a test', 'And so is this']
    db = DataTable(keys, types, descr=descr)
    assert db.data is None, 'No data should be defined'

def test_init():
    keys = ['TEST', 'THIS']
    types = [ int, float]
    descr = ['This is a test', 'And so is this']

    # Test at instantiation
    db = DataTable(keys, types, descr=descr, shape=10)
    assert db.data is not None, 'Data should have been instantiated'
    assert db.shape == (10,), 'Shape is incorrect'

    # Test after the fact
    db = DataTable(keys, types, descr=descr)
    db.init(10)
    assert db.data is not None, 'Data should have been instantiated'
    assert db.shape == (10,), 'Shape is incorrect'

def test_dtype():
    keys = ['TEST', 'THIS']
    types = [ int, float]
    descr = ['This is a test', 'And so is this']
    db = DataTable(keys, types, descr=descr, shape=10)

    _dtype = numpy.dtype([(k,t) for k,t in zip(keys,types)])
    assert db.dtype == _dtype, 'dtype is wrong'

def test_access():
    keys = ['TEST', 'THIS']
    types = [ int, float]
    descr = ['This is a test', 'And so is this']
    db = DataTable(keys, types, descr=descr, shape=10)

    assert numpy.array_equal(db['TEST'], numpy.zeros(10, dtype=int)), \
                'Instantiation array is wrong.'

    db['THIS'] = numpy.arange(10)
    assert numpy.array_equal(db['THIS'], numpy.arange(10, dtype=float)), \
                'Assignment went wrong.'

    db['THIS'][5] = 20.
    assert db['THIS'][5] == 20., 'Assignment of single value went wrong.'

    db[5]['THIS'] = 100.
    assert db['THIS'][5] == 100., 'Assignment of single value went wrong.'

def test_tohdu():
    keys = ['TEST', 'THIS']
    types = [ int, float]
    descr = ['This is a test', 'And so is this']
    db = DataTable(keys, types, descr=descr, shape=10)

    hdu = fits.HDUList([fits.PrimaryHDU(), db.to_hdu(name='TABLE')])

    assert isinstance(hdu['TABLE'], fits.BinTableHDU), 'Bad return type'
    assert numpy.array_equal(hdu['TABLE'].data['TEST'], db['TEST']), 'Arrays are different.'
    assert numpy.array_equal(hdu['TABLE'].data['THIS'], db['THIS']), 'Arrays are different.'

def test_torst():
    keys = ['TEST', 'THIS']
    types = [ int, float]
    descr = ['This is a test', 'And so is this']
    db = DataTable(keys, types, descr=descr, shape=10)
    # Just test for success of the method
    lines = db.to_rst_table()

def test_stellarkin_dt():
    db = StellarKinematicsFitDataTable(ntpl=20, nadd=8, nmult=8, nkin=2, mask_dtype=numpy.int16,
                                       shape=10)
    assert db.shape == (10,), 'Shape is wrong.'
    assert db['TPLWGT'].shape == (10,20), 'Number of template weights is wrong.'
    assert db['ADDCOEF'].shape == (10,8), 'Number of additive polynomial coefficients is wrong'
    assert db['KIN'][0].shape == (2,), 'Number of kinematic moments is wrong.'

def test_emission_dt():
    db = EmissionLineFitDataTable(20, 2, numpy.int16, shape=10)
    assert db.shape == (10,), 'Shape is wrong.'
    assert db.size == 10, 'Size is wrong'
    assert db['BINID'].shape == (10,), 'Number of bin IDs is wrong.'
    assert db['MASK'][0].shape == (20,), 'Number of maskbits is wrong'
    assert db['KIN'].shape == (10,20,2), 'Number of kinematic moments is wrong.'

if __name__ == '__main__':
    test_emission_dt()
