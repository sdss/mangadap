from IPython import embed

import numpy

from mangadap.util.sampling import spectrum_velocity_scale
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.emissionlinetemplates import EmissionLineTemplates
from mangadap.tests.util import data_test_file
from mangadap.contrib.xjmc import ppxf_tied_parameters

def test_init_all():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    dbs = EmissionLineDB.available_databases()
    for key in dbs.keys():
        if 'MSK' in key:
            continue
        print(key)
        emldb = EmissionLineDB.from_key(key)
        etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)


def test_tied():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    emldb = EmissionLineDB(data_test_file('elp_ha_tied.par'))
    etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)

    assert etpl.ntpl == 4, 'Should be 4 templates'
    assert numpy.array_equal(etpl.comp, numpy.arange(etpl.ntpl)), \
            'Each template should be its own kinematic component'
    assert numpy.array_equal(etpl.vgrp, numpy.zeros(etpl.ntpl)), 'All velocities should be tied'
    assert numpy.array_equal(etpl.sgrp, numpy.arange(etpl.ntpl)), 'All dispersions should be untied'


def test_ineq():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    emldb = EmissionLineDB(data_test_file('elp_ha_disp_ineq.par'))
    etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)

    assert etpl.ntpl == 4, 'Should be 4 templates'
    assert numpy.array_equal(etpl.comp, numpy.arange(etpl.ntpl)), \
            'Each template should be its own kinematic component'
    assert numpy.array_equal(etpl.vgrp, numpy.zeros(etpl.ntpl)), 'All velocities should be tied'
    assert numpy.array_equal(etpl.sgrp, numpy.arange(etpl.ntpl)), 'All dispersions should be untied'

    assert etpl.A_ineq is not None, 'Inequality matrix should not be None'
    assert etpl.A_ineq.shape[1] == 2*etpl.ntpl, 'Number of columns in inequality matrix incorrect'


def test_multicomp():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    emldb = EmissionLineDB(data_test_file('elp_ha_multicomp.par'))
    etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)

    assert etpl.ntpl == 12, 'Should be 12 templates'
    assert numpy.amax(etpl.comp)+1 == 10, 'Number of kinematic components changed'
    assert numpy.array_equal(etpl.vgrp, numpy.zeros(etpl.ntpl)), 'All velocities should be tied'
    assert numpy.array_equal(etpl.comp, etpl.sgrp), \
            'Dispersion groups should match kinematic components'

    assert etpl.A_ineq is not None, 'Inequality matrix should not be None'
    assert etpl.A_ineq.shape[1] == 2*(numpy.amax(etpl.comp)+1), \
            'Inequality matrix should have 2 columns per kinematic component'

    ncomp = numpy.amax(etpl.comp)+1
    tied = ppxf_tied_parameters(etpl.comp, etpl.vgrp, etpl.sgrp, numpy.array([2]*ncomp))

    _tied = numpy.array([len(p) > 0 for p in numpy.asarray(tied).flat]).reshape(ncomp,2)

    assert not numpy.any(_tied[:,1]), \
            'No dispersion should be identically tied; all tied line dispersions are either in ' \
            'the same template or same kinematic component.'

    assert numpy.sum(_tied[:,0]) == ncomp-1, 'All but one velocity should be tied.'
    _tied = numpy.asarray(tied)
    assert numpy.all(_tied[1:,0] == _tied[1,0]), \
            'All velocities should be tied to the same component'


def test_multicomp_broadv():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    emldb = EmissionLineDB(data_test_file('elp_ha_multicomp_broadv.par'))
    etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)

    assert etpl.ntpl == 12, 'Should be 12 templates'
    assert numpy.amax(etpl.comp)+1 == 10, 'Number of kinematic components changed'
    assert numpy.array_equal(numpy.unique(etpl.vgrp), [0,1]), 'Should be two velocity groups'
    assert numpy.array_equal(etpl.comp, etpl.sgrp), \
            'Dispersion groups should match kinematic components'

    assert etpl.A_ineq is not None, 'Inequality matrix should not be None'
    assert etpl.A_ineq.shape[1] == 2*(numpy.amax(etpl.comp)+1), \
            'Inequality matrix should have 2 columns per kinematic component'

    ncomp = numpy.amax(etpl.comp)+1
    tied = ppxf_tied_parameters(etpl.comp, etpl.vgrp, etpl.sgrp, numpy.array([2]*ncomp))

    _tied = numpy.array([len(p) > 0 for p in numpy.asarray(tied).flat]).reshape(ncomp,2)

    assert not numpy.any(_tied[:,1]), \
            'No dispersion should be identically tied; all tied line dispersions are either in ' \
            'the same template or same kinematic component.'

    assert numpy.sum(_tied[:,0]) == ncomp-2, 'All but two velocities should be tied.'
    _tied = numpy.asarray(tied)
    assert numpy.all(_tied[2:6,0] == _tied[2,0]), \
            'There should be two groups of tied velocities; this is the first group'
    assert numpy.all(_tied[6:,0] == _tied[6,0]), \
            'There should be two groups of tied velocities; this is the second group'



