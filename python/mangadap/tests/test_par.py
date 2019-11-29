import pytest

from mangadap.par.parset import ParSet

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

#-----------------------------------------------------------------------------

def test_parset():

    # Instantiate a ParSet
    keys = [ 'test', 'par', 'list', 'junk']
    defaults = [ 'this', 0, 0.0, None ]
    options = [ 'this', None, None, None ]
    dtypes = [ str, int, [int, float], None ]
    par = ParSet(keys, defaults=defaults, options=options, dtypes=dtypes)

    # Test it
    assert all([ k in par.keys() for k in keys]), 'Not all keys assigned'
    assert par['junk'] is None, 'Default should be None'
    assert par['test'] == 'this', 'Bad instantiation'

    with pytest.raises(ValueError):
        par['test'] = 'that', 'Options not correctly checked.'
    with pytest.raises(TypeError):
        par['par'] = 'test', 'Type not correctly checked'

    # Try adding a new parameter
    par.add('echo', 10, dtype=int)
    assert 'echo' in par.keys(), 'Key not added'
    with pytest.raises(TypeError):
        par['echo'] = 1.3, 'Type not added'
    
