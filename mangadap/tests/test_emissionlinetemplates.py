from IPython import embed

import numpy

from mangadap.util.sampling import spectrum_velocity_scale
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.emissionlinetemplates import EmissionLineTemplates
from mangadap.tests.util import data_test_file


def test_init_all():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    dbs = EmissionLineDB.available_databases()
    for key in dbs.keys():
        if 'MSK' in key:
            continue
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

