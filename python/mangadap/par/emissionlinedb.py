# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for a database of emission-line parameters, as well as
support classes and functions.

Class usage examples
--------------------

To define an emission line::

    from mangadap.par.emissionlinedb import EmissionLinePar
    p = EmissionLinePar(44, 'Ha', 6564.632, action='f', line='l',
                        flux=1.0, vel=0.0, sig=10., mode='f')

More often, however, you will want to define an emission-line
database using an SDSS parameter file. To do so, you can use one of
the default set of available emission-line databases::

    from mangadap.par.emissionlinedb import EmissionLineDB
    print(EmissionLineDB.available_databases())
    emldb = EmissionLineDB.from_key('ELPMPL9')

The above call uses the :func:`EmissionLineDB.from_key` method to
define the database using its keyword and the database provided with
the MaNGA DAP source distribution. You can also define the database
directly for an SDSS-style parameter file::

    from mangadap.par.emissionlinedb import EmissionLineDB
    emldb = EmissionLineDB('/path/to/emission/line/database/myeml.par')

The above will read the file and set the database keyword to
'MYEML' (i.e., the capitalized root name of the ``*.par`` file).
See :ref:`emissionlines` for the format of the parameter file.

Revision history
----------------

    | **17 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_ and
        `pydl.goddard.astro.airtovac`_ instead of internal functions
    | **13 Jul 2016**: (KBW) Include log_bounded, blueside, and redside
        in database.
    | **06 Oct 2017**: (KBW) Add function to return channel names
    | **02 Dec 2019**: (KBW) Significantly reworked to use the new
        base class.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import os
import glob
import numpy

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny

from .parset import ParSet
from .spectralfeaturedb import SpectralFeatureDB
from ..proc import util

class EmissionLinePar(ParSet):
    r"""
    Parameter object that defines a set of emission-line parameters used
    by various algorithms in the DAP.

    .. todo::
        - Specify these algorithms
        - provide some basic printing functions for user-level
          interaction

    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    Args:
        index (int) : An index used to refer to the line in the *line*
            and *mode* attributes.
        name (str) : A name for the line.
        restwave (float) : The rest wavelength of the line in angstroms
            *in vacuum*.
        action (str) : (**Optional**) Describes how the line should be
            treated; default is ``'f'``.  Possible values are:

            - ``'i'``: ignore the line, as if the line were commented
              out.
            - ``'f'``': fit the line and/or mask the line when fitting
              the stellar continuum.
            - ``'m'``': mask the line when fitting the stellar continuum
              but do NOT fit the line itself
            - ``'s'``: defines a sky line that should be masked.  When
              masked, the wavelength of the line is NOT adjusted for the
              redshift of the object spectrum.

        flux (float) : (**Optional**) **Relative** flux of the emission
            (positive) or absorption (negative) lines.  This should most
            often be unity if ``line='l'`` and indicates the ratio of
            line flux if ``line='dN'``.  Default is 1.0.

        mode (float) : (**Optional**) Fitting mode for the line; default
            is ``'f'``.  Possible values are:

                - ``'f'``: Fit the line independently of all others in
                  its own window.
                - ``'wN'``: Fit the line with untied parameters, but use
                  a window that includes both this line and the line
                  with index N.
                - ``'xN'``: Fit the line with its flux tied to the line
                  with index N.
                - ``'vN'``: Fit the line with the velocity tied to the
                  line with index N.
                - ``'sN'``: Fit the line with the velocity dispersion
                  tied to the line with index N.
                - ``'kN'``: Fit the line with the velocity and velocity
                  dispersion tied to the line with index N.
                - ``'aN'``: Fit the line with the flux, velocity, and
                  velocity dispersion tied to the line with index N.

        profile (str): (**Optional**)  The class definition of the
            profile shape.  This is kept as a string until it is used.
            Once it is used, it is converted to the class name using::
        
                eval(profile)
            
            This line will fail if the profile type has not been
            defined.  Default is ``'GaussianLineProfile'``.
        ncomp (int): (**Optional**) The number of components (number of
            separate line profiles) to use when fitting the line.
            Default is 1.
        output_model (bool): (**Optional**) Flag to include the
            best-fitting model of the line in the emission-line model
            spectrum. Default is True.
        par (numpy.ndarray) : (**Optional**) A list of the initial guess
            for the line profile parameters.  The number of parameters
            must match the struct declaration at the top of the file.
            The initial parameters are automatically adjusted to provide
            any designated flux ratios, and the center is automatically
            adjusted to the provided redshift for the spectrum.  For
            example, for a GaussianLineProfile, this is typically set to
            ``[1.0, 0.0, 100.0]``.
        fix (numpy.ndarray): (**Optional**) A list of flags for fixing
            the input guess parameters during the fit.  Use 0 for a free
            parameter, 1 for a fixed parameter.  The parameter value is
            only fixed **after** adjusted in the flux and or center
            based on the redshift and the implied tied parameters.  For
            a free set of parameters using a GaussianLineProfile, this
            is set to ``[ 0, 0, 0 ]``.
        lobnd (numpy.ndarray): (**Optional**) A list of lower bounds for
            the parameters.  For each parameter, use None to indicate no
            lower bound.  For a GaussianLineProfile with positive flux
            and standard deviation, this is set to ``[ 0.0, None, 0.0
            ]``.
        hibnd (numpy.ndarray): (**Optional**) A list of upper bounds for
            the parameters.  For each parameter, use None to indicate no
            upper bound.  For a GaussianLineProfile with maximum
            standard deviation of 500 km/s, this is set to ``[ None,
            None, 500.0 ]``.
        log_bnd (numpy.ndarray): (**Optional**) A list of flags used
            when determining if a fit parameter is near the imposed
            boundary.  If true, the fraction of the boundary range used
            is done in logarithmic, not linear, separation.
        blueside (numpy.ndarray): (**Optional**) A two-element vector
            with the starting and ending wavelength for a bandpass
            blueward of the emission line, which is used to set the
            continuum level near the emission line when calculating the
            equivalent width.
        redside (numpy.ndarray): (**Optional**) A two-element vector
            with the starting and ending wavelength for a bandpass
            redward of the emission line, which is used to set the
            continuum level near the emission line when calculating the
            equivalent width.
    """
    def __init__(self, index, name, restwave, action=None, flux=None, mode=None, profile=None,
                 ncomp=None, output_model=None, par=None, fix=None, lobnd=None, hibnd=None,
                 log_bnd=None, blueside=None, redside=None):
        
        in_fl = [ int, float ]
        action_options = [ 'i', 'f', 'm', 's']
        arr_like = [ numpy.ndarray, list ]
        
        l = [ name, action, mode ]
        for i in range(len(l)):
            if l[i] is None:
                continue
            l[i] = l[i].strip()

        pars =     [ 'index', 'name', 'restwave', 'action', 'flux', 'mode', 'profile', 'ncomp',
                     'output_model', 'par', 'fix', 'lobnd', 'hibnd', 'log_bnd', 'blueside',
                     'redside' ]
        values =   [   index,   l[0],   restwave,     l[1],   flux,   l[2],   profile,   ncomp,
                       output_model, par,   fix,   lobnd,   hibnd, log_bnd, blueside, redside ]
        defaults = [ None, None, None, 'f', 1.0, 'f', 'GaussianLineProfile', 1, True, None, None,
                     None, None, None, None, None ]
        options =  [ None, None, None, action_options,  None,  None, None, None, None, None, None,
                     None, None, None, None, None ]
        dtypes =   [ int, str, in_fl, str, in_fl, str, str, int, bool, arr_like, arr_like,
                     arr_like, arr_like, arr_like, arr_like, arr_like ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    def _check(self):
        """
        Check the parameter list:

            - Amplitude has to be larger than zero.
            - *mode* must be either ``'f'``, ``'wN'``, ``'xN'``,
              ``'vN'``, ``'sN'``, ``'kN'``, ``'aN'``

        .. todo::
            - Add check to __setitem__()?
            - Add check of profile type

        Raises:
            ValueError: Raised if one of the conditions above are not
                met.
        """
        if self.data['action'] != 'i' and not self.data['flux'] > 0:
            warnings.warn('Emission-line fluxes must be larger than 0.  Ignoring line with ' \
                          'index {0}'.format(self.data['index']))
            self.data['action'] = 'i'
        if self.data['mode'][0] not in ['f', 'w', 'x', 'v', 's', 'k', 'a']:
            raise ValueError('Mode must be either independent (f), fit independently but in the ' \
                             'same window (w), with a tied flux ratio (x), with a tied velocity ' \
                             '(v), with a tied velocity dispersion (s), with both kinematics ' \
                             'tied (k), or with the fluxes and kinematics tied (a).')
        # Check the lengths of all the arrays
        npar = len(self.data['par'])
        if numpy.any(numpy.array([len(self.data['fix']), len(self.data['lobnd']),
                                  len(self.data['hibnd']), len(self.data['log_bnd'])]) != npar):
            raise ValueError('Number of parameters must be the same for par, fix, lobnd, hibnd, '
                             'and log_bnd.')
        if numpy.any(numpy.array([ len(self.data['blueside']), len(self.data['redside'])]) != 2):
            raise ValueError('Bandpasses must be two-element vectors!')


class EmissionLineDB(SpectralFeatureDB):
    r"""
    Basic container class for the database of emission-line parameters.

    See the base class for additional attributes.

    The primary instantiation requires the SDSS parameter file with
    the emission-line data. To instantiate using a keyword (and
    optionally a directory that holds the parameter files), use the
    :func:`mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key`
    class method.

    Args:
        parfile (:obj:`str`):
            The SDSS parameter file with the emission-line database.

    Attributes:
        key (:obj:`str`):
            Database signifying keyword
        file (:obj:`str`):
            File with the emission-line data
        size (:obj:`int`):
            Number of emission lines in the database. 
    """
    default_data_dir = 'emission_lines'
    def _parse_yanny(self):
        """
        Parse the yanny file (provided by :attr:`file`) for the emission-line database.

        Returns:
            :obj:`list`: The list of
            :class:`mangadap.par.parset.ParSet` instances for each
            line of the database.
        """
        # Read the yanny file
        par = yanny(filename=self.file, raw=True)
        if len(par['DAPEML']['index']) == 0:
            raise ValueError('Could not find DAPEML entries in {0}!'.format(self.file))

        # Setup the array of emission line database parameters
        self.size = len(par['DAPEML']['index'])
        parlist = []
        for i in range(self.size):
            invac = par['DAPEML']['waveref'][i] == 'vac'
            # Convert None's to +/- inf
            lobnd = [ -numpy.inf if p == 'None' else float(p) \
                            for p in par['DAPEML']['lower_bound'][i]]
            hibnd = [ numpy.inf if p == 'None' else float(p) \
                            for p in par['DAPEML']['upper_bound'][i] ]
            parlist += [ EmissionLinePar(par['DAPEML']['index'][i], par['DAPEML']['name'][i],
                    par['DAPEML']['lambda'][i] if invac else airtovac(par['DAPEML']['lambda'][i]),
                                         action=par['DAPEML']['action'][i],
                                         flux=par['DAPEML']['relative_flux'][i],
                                         mode=par['DAPEML']['mode'][i],
                                         profile=par['DAPEML']['profile'][i],
                                         ncomp=par['DAPEML']['ncomp'][i],
                                         output_model=bool(par['DAPEML']['output_model'][i]),
                                         par=par['DAPEML']['par'][i],
                                         fix=par['DAPEML']['fix'][i],
                                         lobnd=lobnd, hibnd=hibnd,
                                         log_bnd=par['DAPEML']['log_bounded'][i],
                                blueside=par['DAPEML']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPEML']['blueside'][i])),
                                redside=par['DAPEML']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPEML']['redside'][i])) ) ]
        return parlist

    def channel_names(self, dicttype=True):
        """
        Return a dictionary with the channel names as the dictionary
        key and the channel number as the dictionary value. If
        ``dicttype`` is False, a list is returned with just the
        string keys.
        """
        channels = [ '{0}-{1}'.format(self.data['name'][i], int(self.data['restwave'][i])) 
                            for i in range(self.size) ]
        return { n:i for i,n in enumerate(channels) } if dicttype else channels
