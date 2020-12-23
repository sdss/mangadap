
import pytest
import os

from IPython import embed

import numpy

from mangadap.util.sampling import spectrum_velocity_scale
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.emissionlinetemplates import EmissionLineTemplates


def test_init_all():
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    dbs = EmissionLineDB.available_databases()
    for key in dbs.keys():
        if 'MSK' in key:
            continue
        emldb = EmissionLineDB.from_key(key)
        etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)


