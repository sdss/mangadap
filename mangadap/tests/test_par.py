import pytest

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from astropy.io import fits

from mangadap.par.parset import ParSet
from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.par.bandheadindexdb import BandheadIndexDB
from mangadap.par.absorptionindexdb import AbsorptionIndexDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB

#-----------------------------------------------------------------------------

def test_parset():

    # Instantiate a ParSet
    keys = [ 'test', 'par', 'list', 'junk']
    defaults = [ 'this', 0, 0.0, None ]
    options = [ 'this', None, None, None ]
    dtypes = [ str, int, [int, float], None ]
    descr = [ 'test parameter description',
              'par parameter description',
              None,
              'this is a rather long description of the junk parameter.  You can include ' \
                'rst-style references like pointing back to the ' \
                ':class:`mangadap.par.parset.ParSet` class, for when this description is ' \
                'written to an rst table using :func:`mangadap.par.parset.ParSet.to_rst_table` ' \
                'and included in an rst doc synthesized into html using sphinx.']
    par = ParSet(keys, defaults=defaults, options=options, dtypes=dtypes, descr=descr)

    # Print it
    par.info()
    print(par)

    # Test values
    assert all([ k in par.keys() for k in keys]), 'Not all keys assigned'
    assert par['junk'] is None, 'Default should be None'
    assert par['test'] == 'this', 'Bad instantiation'

    # Test constraints
    with pytest.raises(KeyError):
        par['that'] = 'new', 'This key does not exist; need to use add method.'
    with pytest.raises(ValueError):
        par['test'] = 'that', 'Options not correctly checked.'
    with pytest.raises(TypeError):
        par['par'] = 'test', 'Type not correctly checked'

    # Test IO
    ParSet.from_dict(par.to_dict())

    hdr = fits.Header()
    par.to_header(hdr)
    ParSet.from_header(hdr)

    ParSet.from_config(par.to_config())

    # Try adding a new parameter
    par.add('echo', 10, dtype=int)
    assert 'echo' in par.keys(), 'Key not added'
    with pytest.raises(TypeError):
        par['echo'] = 1.3, 'Type not added'


def test_artifact():

    artdb = ArtifactDB.from_key('BADSKY')
    assert artdb.size == 1

    _artdb = ArtifactDB(artdb.file)
    assert _artdb.size == 1


def test_emldb():
    emldb = EmissionLineDB.from_key('ELPMPL11')
    assert emldb.size == 35

    _emldb = EmissionLineDB(emldb.file)
    assert _emldb.size == 35


def test_eml_channel_names():
    emldb = EmissionLineDB.from_key('ELPMPL11')
    names = emldb.channel_names()
    assert 'Ha-6564' in list(names.keys())
    assert names['OII-3727'] == 0


def test_elmom():
    elmom = EmissionMomentsDB.from_key('ELBMPL9')
    assert elmom.size == 35

    _elmom = EmissionMomentsDB(elmom.file)
    assert _elmom.size == 35


def test_elmom_channel_names():
    elmom = EmissionMomentsDB.from_key('ELBMPL9')
    names = elmom.channel_names()
    assert 'Ha-6564' in list(names.keys())
    assert names['OIId-3728'] == 0


def test_bhddb():

    bhddb = BandheadIndexDB.from_key('BHBASIC')
    assert bhddb.size == 3

    _bhddb = BandheadIndexDB(bhddb.file)
    assert _bhddb.size == 3


def test_bhd_channel_names():
    bhddb = BandheadIndexDB.from_key('BHBASIC')
    names = bhddb.channel_names()
    assert 'D4000' in list(names.keys())
    assert names['D4000'] == 0


def test_absdb():

    absdb = AbsorptionIndexDB.from_key('LICKINDX')
    assert absdb.size == 21

    _absdb = AbsorptionIndexDB(absdb.file)
    assert _absdb.size == 21


def test_abs_channel_names():
    absdb = AbsorptionIndexDB.from_key('LICKINDX')
    names = absdb.channel_names()
    assert 'CN1' in list(names.keys())
    assert names['CN1'] == 0





