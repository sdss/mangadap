
.. include:: include/links.rst

.. _gettingstarted:

Getting started
===============

This is your "Quick-Start" guide to working with the DAP output and
serves as an introduction to the other, more detailed pages
discussing the DAP output. Before diving into the MaNGA data, you are
strongly encouraged to briefly familiarize yourself with these pages,
in particular.

 * :ref:`metadatamodel`
 * :ref:`datamodel`
 * :ref:`qualityassessment`
 * :ref:`corrections`
 * :ref:`aggregation`

The DAP provides:

 * Spaxel-by-spaxel coordinate information based on the isophotal
   ellipticities from the NASA-Sloan Atlas elliptical Petrosian analysis
 * S/N measurements
 * Binned spectra
 * Stellar kinematics and stellar-continuum models
 * Emission-line properties and models
 * Spectral-index measurements

It is important that you understand the procedures and algorithms
used by the DAP to provide the data you're interested in if you are
to trust their scientific usage. In particular, consult the DAP
Overview paper, `Westfall et al. (2019, AJ, 158, 231)`_, the
Emission-Line Modeling paper, `Belfiore et al. (2019, AJ, 158,
160)`_, and the :ref:`workflow`.

.. note::

    This page is meant to provide a useful guide for how to get
    started with the DAP output data. *If anything is unclear to you
    then this page is not doing its job!*

    Please `Submit an issue`_ if you have any questions.

    Also, if you have helpful code snippets that you'd like to share,
    please do! The best way to do so would be via a pull request.
    Please include your code as a script in the ``examples``
    directory, and include any necessary modules in the
    ``mangadap/contrib`` directory. If necessary, please include
    licenses relevant to your code in your PR.

----

Directory structure
-------------------

The :ref:`datamodel-directory-structure` is described as part of the
full DAP :ref:`datamodel`. In particular, note the root directories to
the:

 * :ref:`metadatamodel-dapall`,
 * :ref:`datamodel-maps`, and
 * :ref:`datamodel-cube`

----

.. _gettingstarted-daptype:

DAPTYPE selection
-----------------

The DAP analyzes each DRP-produced datacube multiple ways. Each
:ref:`datamodel-daptype` is given a unique keyword, called the
``DAPTYPE``. 

The ``DAPTYPEs`` in DR15 are:

 * ``VOR10-GAU-MILESHC``: Analysis of spectra binned to :math:`{\rm
   S/N}\sim10` using the Voronoi binning algorithm `Cappellari & Copin (2003,
   MNRAS, 342, 345)`_

 * ``HYB10-GAU-MILESHC``: Stellar-continuum analysis of spectra
   binned to :math:`{\rm S/N}\sim10` for the stellar kinematics (same as
   the ``VOR10`` approach); however, the emission-line measurements are
   performed on the individual spaxels.

The ``DAPTYPEs`` in MPL-11 are:

 * ``SPX-MILESHC-MASTARSSP``: Analysis of each individual spaxel; spaxels
   must have a valid continuum fit for an emission-line model to be fit

 * ``VOR10-MILESHC-MASTARSSP``: Analysis of spectra binned to :math:`{\rm
   S/N}\sim10` using the Voronoi binning algorithm `Cappellari & Copin (2003,
   MNRAS, 342, 345)`_

 * ``HYB10-MILESHC-MASTARSSP``: Stellar-continuum analysis of spectra
   binned to :math:`{\rm S/N}\sim10` for the stellar kinematics (same
   as the ``VOR10`` approach); however, the emission-line and
   spectral-index measurements are performed on the individual
   spaxels.

 * ``HYB10-MILESHC-MASTARHC2``: Same as the above except the hierarchically
   clustered MaStar stellar spectra are used to fit the stellar continuum in
   the emission-line module.

The type of analysis that you should use will depend on your science
application. In all cases, please consult `Westfall et al. (2019, AJ,
158, 231)`_ and `Belfiore et al. (2019, AJ, 158, 160)`_ for usage
guidelines and limitations of the data.

SPX-MILESHC-MASTARSSP (MPL-11 only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are useful for most science applications that can push to very
low S/N. They are also useful for characterizing the performance of
the measurements toward the low S/N limit of the data.

Important considerations:

 * Spectra with :math:`g`-band :math:`{\rm S/N} < 1` will not have a
   stellar-continuum model or Gaussian emission-line models.

 * Beware of biases in the stellar velocity dispersion with :math:`{\rm
   S/N}_g < 10`; cf. `Westfall et al. (2019, AJ, 158, 231)`_

VOR10-GAU-MILESHC (DR15), VOR10-MILESHC-MASTARSSP (MPL-11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These data are geared toward stellar kinematics, where the data are
Voronoi binned to a :math:`g`-band :math:`{\rm S/N}\sim 10`.

Important considerations:

 * Spectra with :math:`g`-band :math:`{\rm S/N} < 1` are not included
   in any bin.
 * Voronoi-binned spectra are just simple means of all the spectra in
   the bin.
 * The covariance in the datacube is propagated to the variance in
   the stacked spectra.
 * The spectral resolution in each binned spectra is propagated to
   the expected spectral resolution in the stacked spectrum, similar
   to the formalism used by
   :func:`~mangadap.spectra.rowstackedspectra.RowStackedSpectra.instrumental_dispersion_plane`.
 * (Binned) Spectra with :math:`g`-band :math:`{\rm S/N} < 1` will
   *not* have a stellar-continuum model or Gaussian emission-line
   model.
 * Because the binning is done based on the *continuum* S/N, this
   limits the emission-line science that can be done at low continuum
   S/N.

HYB10-GAU-MILESHC (DR15), HYB10-MILESHC-MASTARSSP (MPL-11), HYB10-MILESHC-MASTARHC2 (MPL-11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the default files that most users will want to use. We
first Voronoi-binned the spectra to a :math:`g`-band :math:`{\rm
S/N}\sim 10` to measure the stellar kinematics. Then these bins are
deconstructed to fit the emission lines. We refer to this approach as
the :ref:`datamodel-hybrid-binning`.

Important considerations:

 * Spectra with :math:`g`-band :math:`{\rm S/N} < 1` are not included
   in any bin.
 * Voronoi-binned spectra are just simple means of all the spectra in
   the bin.
 * The covariance in the datacube is propagated to the variance in
   the stacked spectra.
 * The spectral resolution in each binned spectra is propagated to
   the expected spectral resolution in the stacked spectrum, similar
   to the formalism used by
   :func:`~mangadap.spectra.rowstackedspectra.RowStackedSpectra.instrumental_dispersion_plane`.
 * (Binned) Spectra with :math:`g`-band :math:`{\rm S/N} < 1` will
   not have a stellar-continuum model.
 * All spectra with 80% valid pixels will have a combined
   emission-line+stellar-continuum model, where the stellar
   kinematics have been fixed by the fits to the (nearest) binned
   spectrum.
 * This is the only file where the bin IDs are different for the
   emission-line properties and spectral indices.

----

Output files
------------

The primary output files are located at:

+---------------+
| SAS Directory |
+===============+
|   `DR15 SAS`_ |
+---------------+
| `MPL-11 SAS`_ |
+---------------+

There are two main output files for each observation
(``plate``-``ifudesign`` combination):

 * ``manga-[PLATE]-[IFUDESIGN]-MAPS-[DAPTYPE].fits.gz``, see
   :ref:`datamodel-maps`: 2D "maps" (i.e., images) of DAP measured
   properties.
 * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE].fits.gz``, see
   :ref:`datamodel-cube`: 3D data cubes with the binned and
   best-fitting-model spectra.

The datacubes produced by the DAP have the same shape as the DRP
datacube, and the DAP maps have the same spatial dimensions as a single
wavelength channel in the DRP datacubes.  This is meant to ease
associating the DRP input and DAP output products.

Examples are given below for how to interact with the two main output
files using python. However, you are **strongly encouraged** to
`install Marvin`_ and use it to interact with the data.

----

Additionally, the DAP produces a summary file, the
:ref:`metadatamodel-dapall`, with global metrics for each observation
and analysis approach (DAPTYPE). This file is at the top-level
directory of each MPL/DR.

----

MAPS files
~~~~~~~~~~

:ref:`datamodel-maps`: The ``MAPS`` files are the primary output file
from the DAP.

In brief, the file contains 2D "maps" (i.e., images) of DAP measured
properties.  Most properties are provided in groups of three fits
extensions:

  #. ``[property]``: the measurement value,
  #. ``[property]_IVAR``: the measurement uncertainty stored as the
     inverse variance, and
  #. ``[property]_MASK``: a corresponding bit mask for each spaxel.

The headers of each extension provides the astrometric World Coordinate
System (WCS) and should exactly match that of the DRP output
``LOGCUBE`` files (apart from the wavelength coordinate).

Many properties have multiple "species" or channels associated with
them. The identifying name of each mapped property is provided in the
header; e.g., the emission lines fit in MPL-11 are listed :ref:`here
<datamodel-emission-line-channels>` as they appear in the header of
the relevant ``MAPS`` file extension (e.g., ``EMLINE_GFLUX``). In
python, you can create a dictionary of items in each channel using
the following method:

.. code-block:: python 

    # Declare a function that creates a dictionary for the columns in the
    # multi-channel extensions
    def channel_dictionary(hdu, ext, prefix='C'):
        """
        Construct a dictionary of the channels in a MAPS file.
        """
        channel_dict = {}
        for k, v in hdu[ext].header.items():
            if k[:len(prefix)] == prefix:
                try:
                    i = int(k[len(prefix):])-1
                except ValueError:
                    continue
                channel_dict[v] = i
        return channel_dict

    def channel_units(hdu, ext, prefix='U'):
        """
        Construct an array with the channel units.
        """
        cu = {}
        for k, v in hdu[ext].header.items():
            if k[:len(prefix)] == prefix:
                try:
                    i = int(k[len(prefix):])-1
                except ValueError:
                    continue
                cu[i] = v.strip()
        channel_units = numpy.empty(max(cu.keys())+1, dtype=object)
        for k, v in cu.items():
            channel_units[k] = v
        return channel_units.astype(str)


This is identical to :func:`~mangadap.util.fileio.channel_dictionary`,
such that you can do the following:

.. code-block:: python

    from mangadap.util.fileio import channel_dictionary
    from astropy.io import fits

    hdu = fits.open('manga-8138-12704-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz')
    emlc = channel_dictionary(hdu, 'EMLINE_GFLUX')
    ha_flux = hdu['EMLINE_GFLUX'].data[emlc['Ha-6564']]

In general, it is best to **select the extension and channel based on
its name**, *not* its number. This is because the ordering of, e.g.,
the emission lines in the relevant extensions has changed between
different DRs/MPLs and may change again. An example of how to do this
is shown :ref:`below <gettingstarted-maps-example>`.

.. _gettingstarted-maps-example:

Usage example
+++++++++++++

With the ``MAPS`` fits file, you should be able to extract DAP maps
output using any fits reader. **Please** `Submit an issue`_ if you run
into any problems!

For example, we provide a python code snippet below that will plot the
:math:`{\rm H}\alpha` flux map, stellar velocity field, the corrected
stellar velocity dispersion field, and the corrected :math:`{\rm
H}\beta` index map for
`manga-8138-12704-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz`_.
Importantly, note the use of the :ref:`corrections`!!

.. code-block:: python

    # Imports
    import numpy
    from matplotlib import pyplot
    from astropy.io import fits
    from mangadap.util.fileio import channel_dictionary, channel_units

    def apply_index_dispersion_correction(indx, indxcorr, unit):
        """
        Apply a set of dispersion corrections.
        """
        if unit not in [ 'ang', 'mag', '' ]:
            raise ValueError('Unit must be mag or ang.')
        return indx * indxcorr if unit in ['ang',''] else indx + indxcorr

    # Open the fits file
    hdu = fits.open('manga-8138-12704-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz')

    # Build dictionaries with the emission-line and spectral-index
    # channel names to ease selection and get the spectral-index units
    emlc = channel_dictionary(hdu, 'EMLINE_GFLUX')
    spic = channel_dictionary(hdu, 'SPECINDEX')
    spiu = channel_units(hdu, 'SPECINDEX')

    # Show the Gaussian-fitted H-alpha flux map
    mask_ext = hdu['EMLINE_GFLUX'].header['QUALDATA']
    halpha_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['Ha-6564']],
                                       mask=hdu[mask_ext].data[emlc['Ha-6564']] > 0)

    pyplot.imshow(halpha_flux, origin='lower', interpolation='nearest', cmap='inferno')
    pyplot.colorbar()
    pyplot.show()

    # Show the stellar velocity field
    mask_ext = hdu['STELLAR_VEL'].header['QUALDATA']
    stellar_vfield = numpy.ma.MaskedArray(hdu['STELLAR_VEL'].data, mask=hdu[mask_ext].data > 0)

    pyplot.imshow(stellar_vfield, origin='lower', interpolation='nearest', vmin=-300, vmax=300,
                  cmap='RdBu_r')
    pyplot.colorbar()
    pyplot.show()

    # Show the corrected stellar velocity dispersion field
    mask_ext = hdu['STELLAR_SIGMA'].header['QUALDATA']
    stellar_sfield_sqr = numpy.ma.MaskedArray(numpy.square(hdu['STELLAR_SIGMA'].data)
                                              - numpy.square(hdu['STELLAR_SIGMACORR'].data[0]),
                                              mask=hdu[mask_ext].data > 0)
    # WARNING: This will ignore any data where the correction is larger than the measurement
    stellar_sfield = numpy.ma.sqrt(stellar_sfield_sqr)

    pyplot.imshow(stellar_sfield, origin='lower', interpolation='nearest', cmap='viridis')
    pyplot.colorbar()
    pyplot.show()

    # Show the corrected H-beta index map
    mask_ext = hdu['SPECINDEX'].header['QUALDATA']
    hbeta_index_raw = numpy.ma.MaskedArray(hdu['SPECINDEX'].data[spic['Hb']],
                                           mask=hdu[mask_ext].data[spic['Hb']] > 0)
    hbeta_index = apply_index_dispersion_correction(hbeta_index_raw,
                                                    hdu['SPECINDEX_CORR'].data[spic['Hb']],
                                                    spiu[spic['Hb']])
    
    pyplot.imshow(hbeta_index, origin='lower', interpolation='nearest', cmap='inferno')
    pyplot.colorbar()
    pyplot.show()

----

Model LOGCUBE files
~~~~~~~~~~~~~~~~~~~

:ref:`datamodel-cube`: The ``LOGCUBE`` files provide the binned
spectra and the best-fitting model spectrum for each spectrum that was
successfully fit.

These files are useful for detailed assessments of the model
parameters because they allow you to return to the spectra and
compare the model against the data. As described by `Westfall et al.
(2019, AJ, 158, 231)`_, the DAP fits the spectra in two stages, one
to get the stellar kinematics and the second to determine the
emission-line properties. The emission-line module (used for all
binning schemes) fits both the stellar continuum and the emission
lines at the same time, where the stellar kinematics are fixed by the
first fit. The stellar-continuum models from the first fit are
provided in the ``STELLAR`` extension; to get the stellar continuum
determined during the emission-line modeling, you have to subtract
the emission-line model (in the ``EMLINE`` extension) from the full
model (in the ``MODEL`` extension); see our
:ref:`gettingstarted-cube-example`.

.. warning::

    In the ``HYB`` binning case the binned spectra provided in the
    ``LOGCUBE`` files are from the Voronoi binning step. However, the
    emission-line models are fit to the *individual spaxels*. So:

        - The stellar-continuum fits from the first iteration, in the
          ``STELLAR`` extension, should be compared to the Voronoi
          binned spectra in the file, but
        - the best-fitting model spectra in the ``MODEL`` extension
          should be compared to the individual spectra from the DRP
          ``LOGCUBE`` file!

.. _gettingstarted-cube-example:

Usage example
+++++++++++++

With the ``LOGCUBE`` fits file, you should be able to extract the
binned spectra and best-fitting models produced by the DAP using any
fits reading software. **Please** `Submit an issue`_ if you run into
any problems!

For example, we provide a python code snippet below that plots the
highest S/N spectrum, the full model, the residuals, the model
stellar continuum, and the model emission-line spectrum using
`manga-8138-12704-LOGCUBE-HYB10-MILESHC-MASTARSSP.fits.gz`_.
Importantly, note the use of the :ref:`metadatamodel-maskbits`!!

.. code-block:: python

    # Imports
    import numpy
    from astropy.io import fits
    from matplotlib import pyplot

    # This is a bitmask handling object from the DAP source code
    from mangadap.dapfits import DAPCubeBitMask

    # Open the fits file
    hdu_maps = fits.open('manga-8138-12704-MAPS-SPX-MILESHC-MASTARSSP.fits.gz')
    hdu_cube = fits.open('manga-8138-12704-LOGCUBE-SPX-MILESHC-MASTARSSP.fits.gz')

    # Get the S/N per bin from the MAPS file
    snr = numpy.ma.MaskedArray(hdu_maps['BIN_SNR'].data, mask=hdu_maps['BINID'].data[0] < 0)

    # Select the bin/spaxel with the highest S/N and get the relevant map coordiantes
    k = numpy.ma.argmax(snr.ravel())
    j, i = numpy.unravel_index(k, hdu_maps['BIN_SNR'].data.shape)

    # Declare the bitmask object to mask selected pixels
    bm = DAPCubeBitMask()
    wave = hdu_cube['WAVE'].data
    flux = numpy.ma.MaskedArray(hdu_cube['FLUX'].data[:,j,i],
                                mask=bm.flagged(hdu_cube['MASK'].data[:,j,i],
				                                ['IGNORED', 'FLUXINVALID', 'IVARINVALID',
                                                 'ARTIFACT']))

    model = numpy.ma.MaskedArray(hdu_cube['MODEL'].data[:,j,i],
                                 mask=bm.flagged(hdu_cube['MODEL_MASK'].data[:,j,i], 'FITIGNORED'))
    stellarcontinuum = numpy.ma.MaskedArray(hdu_cube['MODEL'].data[:,j,i]
                                                - hdu_cube['EMLINE'].data[:,j,i],
                                            mask=bm.flagged(hdu_cube['MODEL_MASK'].data[:,j,i],
                                                            'FITIGNORED'))
    emlines = numpy.ma.MaskedArray(hdu_cube['EMLINE'].data[:,j,i],
                                   mask=bm.flagged(hdu_cube['MODEL_MASK'].data[:,j,i], 'ELIGNORED'))
    resid = flux-model-0.5

    pyplot.step(wave, flux, where='mid', color='k', lw=0.5)
    pyplot.plot(wave, model, color='r', lw=1)
    pyplot.plot(wave, stellarcontinuum, color='g', lw=1)
    pyplot.plot(wave, emlines, color='b', lw=1)
    pyplot.step(wave, resid, where='mid', color='0.5', lw=0.5)
    pyplot.show()

----

DAPall file
~~~~~~~~~~~

:ref:`metadatamodel-dapall`: The DAP provides a single summary file
for all analyses done for a given MPL/DR.

The name of the file is ``dapall-$MANGADRP_VER-$MANGADAP_VER.fits``;
for example, this is ``dapall-v3_1_1-3.1.0.fits`` in MPL-11.

This file is used primarily for sample selection, and it has many
limitations (e.g., it provides star-formation rates calculated
*without* applying an internal extinction correction).

The code below provides an example use, plotting the stellar mass
from the `DRPall`_ file against the star-formation rates from the
DAPall file.

.. code-block:: python

    from matplotlib import pyplot
    from astropy.io import fits

    drpall = fits.open('drpall-v3_1_1.fits')['MANGA'].data
    daptype = 'HYB10-MILESHC-MASTARSSP'
    dapall = fits.open('dapall-v3_1_1-3.1.0.fits')[daptype]
    mass = drpall['nsa_sersic_mass'][dapall['DRPALLINDX']]

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPQUAL'] == 0) & (mass > 0) \
                & (dapall['SFR_TOT'] > -999)

    pyplot.scatter(mass[indx], dapall['SFR_TOT'][indx], marker='.', s=20, lw=0, color='k')
    pyplot.xlim([1e8, 1e12])
    pyplot.ylim([1e-4, 1e1])
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel(r'$\mathcal{M}_{\ast}$ [$h^{-2} \mathcal{M}_{\odot}$]')
    pyplot.ylabel(r'SFR [$h^{-2} \mathcal{M}_{\odot}$/yr]')
    pyplot.show()

----

.. _gettingstarted-bitmasks:

Using the pixel/spaxel masks
----------------------------

The maskbits for the DAP data are described :ref:`here
<metadatamodel-maskbits>`. In particular, be aware of the
``DONOTUSE`` and the ``UNRELIABLE`` flags for the MAPS files. The 2d
``MAPS`` file pixel mask is :ref:`metadatamodel-dappixmask`. The 3d
``LOGCUBE`` file spaxel mask is :ref:`metadatamodel-dapspecmask`.

In all cases, the DAP has a convenience class that allows a user to
quickly determine if any mask value is flagged with a certain value.
For example:

.. code-block:: python

    # Imports
    import os
    from astropy.io import fits
    from mangadap.util.bitmask import BitMask
    from mangadap.config.defaults import sdss_maskbits_file

    # Define the path to the IDLUTILS maskbits file
    sdssMaskbits = sdss_maskbits_file()

    # Instantiate the BitMask object
    bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_DAPQUAL')

    # Read a DAP file
    hdu = fits.open('manga-8138-12704-MAPS-SPX-MILESHC-MASTARSSP.fits.gz')

    # Check if the file is critical and print the result
    dap_file_is_critical = bm.flagged(hdu['PRIMARY'].header['DAPQUAL'], flag='CRITICAL')
    print('This DAP file {0} flagged as CRITICAL.'.format('is' if dap_file_is_critical
                                                          else 'is not'))

There are also a predefined set of derived
:class:`~mangadap.util.bitmask.BitMask` classes that the DAP provides.
For example:

.. code-block:: python

    #Imports
    import numpy
    from astropy.io import fits
    from mangadap.util.drpfits import DRPFitsBitMask

    # Instantiate the BitMask object
    bm = DRPFitsBitMask()

    # Read a DRP file
    hdu = fits.open('manga-8138-12704-LOGCUBE.fits.gz')

    # Find the number of pixels flagged as DONOTUSE or FORESTAR
    indx = bm.flagged(hdu['MASK'].data, flag=['DONOTUSE', 'FORESTAR']) 
    print('This DRP file has {0}/{1} pixels flagged as either DONOTUSE or FORESTAR.'.format(
            numpy.sum(indx), numpy.prod(indx.shape)))

See also the `Marvin Maskbits utility`_.

----

.. _gettingstarted-binid:

Using the BINID extension
-------------------------

The ``BINID`` extension has 5 channels.  They provide the IDs of spaxels
associated with:

    0. each binned spectrum.  Any spaxel with ``BINID=-1`` as not included in any bin.
    1. any binned spectrum with an attempted stellar kinematics fit.
    2. any binned spectrum with emission-line moment measurements.
    3. any binned spectrum with an attempted emission-line fit.
    4. any binned spectrum with spectral-index measurements.

In any of these channels, you can obtain the unique bin numbers using
``numpy.unique(bin_indx.ravel())[1:]``; the selection of all but the
first array element is just providing all the numbers without the -1 
for invalid spaxels (assuming all bin ID maps will have spaxels that
are not within the IFU field-of-view, which is always true for
MaNGA). If you're working with anything but the ``SPX`` binning,
you'll want to extract the unique spectra and/or maps values. You can
do that by finding the indices of the unique bins, like this:

.. code-block:: python

    unique_bins, unique_indices = tuple(map(lambda x : x[1:], numpy.unique(bin_indx.ravel(),
                                                                           return_index=True)))

Here's a worked example where I use
:func:`~mangadap.util.fitsutil.DAPFitsUtil.unique_bins` to pull out the
unique stellar velocities and produce a scatter plot of the x and y
positions of the luminosity-weighted bin centers and color them by the
measure stellar velocity.

.. code-block:: python

    #Imports
    import numpy
    from astropy.io import fits
    from matplotlib import pyplot

    from mangadap.util.fitsutil import DAPFitsUtil

    # Read a DAP MAPS file
    hdu = fits.open('manga-8138-12704-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz')

    # Get the unique indices of the stellar-kinematics bins
    ubins, uindx = DAPFitsUtil.unique_bins(hdu['BINID'].data[1], return_index=True)

    # Get the x and y coordinates and the stellar velocities
    x = hdu['BIN_LWSKYCOO'].data[0].ravel()[uindx]
    y = hdu['BIN_LWSKYCOO'].data[1].ravel()[uindx]
    v = numpy.ma.MaskedArray(hdu['STELLAR_VEL'].data.ravel()[uindx],
                             mask=hdu['STELLAR_VEL_MASK'].data.ravel()[uindx] > 0)

    fig = pyplot.figure(figsize=pyplot.figaspect(1))

    ax = fig.add_axes([0.15, 0.15, 0.65, 0.65], facecolor='0.8')
    cax = fig.add_axes([0.81, 0.15, 0.02, 0.65])
    ax.minorticks_on()
    ax.set_xlim([18,-18])
    ax.set_ylim([-18,18])
    ax.grid(True, which='major', color='0.7', zorder=0, linestyle='-')

    sp = ax.scatter(x, y, c=v, vmin=-300, vmax=300, cmap='RdBu_r', marker='.', s=30, lw=0, zorder=3)
    pyplot.colorbar(sp, cax=cax)

    ax.text(0.5, -0.1, r'$\xi$ (arcsec)',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.13, 0.5, r'$\eta$ (arcsec)', rotation='vertical',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    cax.text(5, 0.5, r'$V_\ast$ (km/s)', rotation='vertical',
             horizontalalignment='center', verticalalignment='center', transform=cax.transAxes)

    pyplot.show()


