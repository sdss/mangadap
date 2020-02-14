
import pytest
import numpy

from mangadap.util.extinction import GalacticExtinction

def test_ccm():
    ccm = GalacticExtinction(form='CCM')

    Av = 0.068
    ebv = Av/ccm.rv

    wave = numpy.array([ 3500, 5000, 7000, 10000])

    assert numpy.allclose(ccm.compute(wave, ebv),
                          numpy.array([1.10482203, 1.07281569, 1.04810482, 1.02562548])), \
                'Bad CCM redenning values.'

def test_odonnell():
    ccm = GalacticExtinction(form='CCM')
    odonnell = GalacticExtinction(form='ODonnell')

    Av = 0.068
    ebv = Av/odonnell.rv

    wave = numpy.array([ 3500, 5000, 7000, 10000])

    assert numpy.all(odonnell.compute(wave, ebv)/ccm.compute(wave, ebv) < 1.001), \
            'O\'Donnell values incorrect.'

def test_fitzpatrick():
    ccm = GalacticExtinction(form='CCM')
    fitzpatrick = GalacticExtinction(form='FM')

    Av = 0.068
    ebv = Av/fitzpatrick.rv

    wave = numpy.array([ 3500, 5000, 7000, 10000])

    assert numpy.all(fitzpatrick.compute(wave, ebv)/ccm.compute(wave, ebv) < 1.), \
            'Fitzpatrick values incorrect.'

def test_calzetti():
    ccm = GalacticExtinction(form='CCM')
    calzetti = GalacticExtinction(form='Calzetti')

    Av = 0.068
    ebv = Av/calzetti.rv

    wave = numpy.array([ 3500, 5000, 7000, 10000])

    assert numpy.all(calzetti.compute(wave, ebv)/ccm.compute(wave, ebv) < 1.02), \
            'Calzetti values incorrect.'

