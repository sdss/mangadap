# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Defines a catch-all class for useful constants.  These are meant to be
values that are **not** available elsewhere, such as
`astropy.constants`_.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy

class DAPConstants:
    r"""
    Defines the following set of constants:

    +-----------------------------+-------------------------------+
    | Attribute                   | Value                         |
    +=============================+===============================+
    | :attr:`sig2fwhm`            | :math:`\sqrt{8\ln(2)}`        |
    +-----------------------------+-------------------------------+
    | :attr:`rad2arcs`            | :math:`3600\frac{180}{\pi}`   |
    +-----------------------------+-------------------------------+
    | :attr:`sidereal_year`       | :math:`31558175.779`          |
    +-----------------------------+-------------------------------+
    """

    # Convert from sigma to FWHM: FWHM = sig2fwhm * sig
    sig2fwhm = numpy.sqrt(8.0 * numpy.log(2.0))
    r"""
    Conversion factor from the standard deviation, :math:`\sigma`, of
    a Gaussian to its full-width at half maximum (FWHM).
    """

    # Convert from radians to arcseconds: arcsec = rad2arcs * radians
    rad2arcs = 60*60*180/numpy.pi
    """
    Conversion factor from radians to to arcseconds.
    """


    # Length of one sidereal year in seconds ()
    sidereal_year = 31558175.779
    """
    Length of a sidereal year (1.0000385 Gregorian years) in seconds.
    """

    def __init__(self):
        self.sig2fwhm = constants.sig2fwhm
        self.rad2arcs = constants.rad2arcs
        self.sidereal_year = constants.sidereal_year


