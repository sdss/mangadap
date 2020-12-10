
import pytest
import os

from IPython import embed

import numpy

from mangadap.util.sampling import spectrum_velocity_scale
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.par.emissionlinedb import EmissionLineDBNew
from mangadap.proc.emissionlinetemplates import EmissionLineTemplates
from mangadap.proc.emissionlinetemplates import EmissionLineTemplatesNew

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)


def test_init_mpl9():
    emldb = EmissionLineDB.from_key('ELPMPL9')
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    # Set the dispersion to be 1 pixel
    velscale = spectrum_velocity_scale(wave)
    etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)


def test_init_mpl11():
    emldb = EmissionLineDBNew.from_key('ELPMPL11')
    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    # Set the dispersion to be 1 pixel
    velscale = spectrum_velocity_scale(wave)
    etpl = EmissionLineTemplatesNew(wave, velscale, emldb=emldb)

    embed()
    exit()


if __name__ == '__main__':

    test_init_mpl11()

    wave = numpy.logspace(*numpy.log10([3600., 10000.]), 4563)
    velscale = spectrum_velocity_scale(wave)

    old_emldb = EmissionLineDB.from_key('ELPMILES')
    old_etpl = EmissionLineTemplates(wave, velscale, emldb=old_emldb)

    emldb = EmissionLineDBNew.from_key('ELPMILES_NEW')
    etpl = EmissionLineTemplatesNew(wave, velscale, emldb=emldb)

    embed()
    exit()

#    test_init_mpl9()
#    test_init_mpl11()


