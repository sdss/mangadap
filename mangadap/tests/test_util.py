
from IPython import embed
import numpy as np
import pytest
from mangadap.proc import util

def test_covariance_calibrated_noise():

    nbin = 10
    par = 1.62

    noise = util.covariance_calibrated_noise(nbin, par)
    assert noise == (1 + par), 'Bad calculation'

    par = [1.1, 1.62, 100]
    noise = util.covariance_calibrated_noise(nbin, par)
    assert noise == par[1]*(1 + par[0]), 'Bad calculation'

    # Should fail if there are only two parameters
    par = [1.1, 1.62]
    with pytest.raises(ValueError):
        noise = util.covariance_calibrated_noise(nbin, par)

    # Test the influence of the threshold
    nbin = [10, 75]
    par = [1.1, 1.62, 100]
    noise = util.covariance_calibrated_noise(nbin, par)
    assert isinstance(noise, np.ndarray), \
                'Should be an array on output because multiple nbin were passed.'
    
    par = [1.1, 1.62, 50]
    _noise = util.covariance_calibrated_noise(nbin, par)
    assert noise[0] == _noise[0], 'calibrated noise below nbin threshold should be the same'
    assert noise[1] > _noise[1], 'calibrated noise above nbin threshold should be lowered'
    assert _noise[1] == util.covariance_calibrated_noise(par[2], par), 'should be the same'

