
import pytest
import os

from IPython import embed

import numpy

from mangadap.par.absorptionindexdb import AbsorptionIndexDB
from mangadap.par.bandheadindexdb import BandheadIndexDB
from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.spectralindices import SpectralIndices, SpectralIndicesBitMask

import warnings
warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

def test_model_indices():

    # Setup
    bm = SpectralIndicesBitMask()
    absdb = AbsorptionIndexDB.from_key('EXTINDX')
    bhddb = BandheadIndexDB.from_key('BHBASIC')

    # Grab the model spectra
    tpl = TemplateLibrary('M11MILES', spectral_step=1e-4, log=True, hardcopy=False)
    flux = numpy.ma.MaskedArray(tpl['FLUX'].data,
                                mask=tpl.bitmask.flagged(tpl['MASK'].data,
                                                         flag=['NO_DATA', 'WAVE_INVALID',
                                                               'FLUX_INVALID', 'SPECRES_NOFLUX']))

    # Try to measure only absorption-line indices
    indices = SpectralIndices.measure_indices(absdb, None, tpl['WAVE'].data, flux[:2,:],
                                              bitmask=bm)

    # Try to measure only bandhead indices
    indices = SpectralIndices.measure_indices(None, bhddb, tpl['WAVE'].data, flux[:2,:],
                                              bitmask=bm)

    # Measure both
    indices = SpectralIndices.measure_indices(absdb, bhddb, tpl['WAVE'].data, flux[:2,:],
                                              bitmask=bm)

    # Test the output
    indx = numpy.ma.MaskedArray(indices['INDX'], mask=bm.flagged(indices['MASK']))
    indx_bf = numpy.ma.MaskedArray(indices['INDX_BF'], mask=bm.flagged(indices['MASK']))

    assert indx.shape == (2,46), 'Incorrect output shape'
    assert numpy.sum(indx.mask[0,:]) == 12, 'Incorrect number of masked indices'
    assert numpy.allclose(indx[1,:].compressed(),
                numpy.array([-2.67817275e-02,   1.94113040e-02,   2.92495672e-01,
                              1.91132126e+00,   1.86182351e+00,   1.13929718e+00,
                              2.45269081e+00,   2.24190385e+00,   3.55293194e+00,
                              5.25612762e+00,   2.09506842e-02,   5.47343625e-02,
                              2.29779119e-01,   1.92039303e+00,   1.95663049e+00,
                              9.81723064e-01,   7.44058546e-01,   5.46484409e-01,
                              1.56520293e+00,   8.93716574e-03,   1.82158534e-02,
                              2.46748795e+00,   1.03916389e+00,   2.27668045e+00,
                              2.68943564e+00,   6.77715967e+00,   1.31272416e+00,
                              9.12401084e-03,   3.52657124e-03,   1.92826953e-03,
                             -1.00946682e-02,   2.01732938e-02,   1.17072128e+00,
                              1.12525642e+00]),
                          rtol=0.0, atol=1e-4), 'Index values are different'
    assert numpy.std(indx.compressed() - indx_bf.compressed()) < 0.01, \
            'Index definitions are too different'
