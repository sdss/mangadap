

# 10. MaNGA TRM: Getting started with DAP output data #

'''This is your Quick Start guide to working with the DAP products.'''  However, you are strongly encouraged to familiarize yourself with the rest of the DAP documentation in the TRM, particularly with regard to the '''[wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel Data Model]'''.

The DAP provides:
 - Spaxel-by-spaxel coordinate information based on the isophotal ellipticities from the NASA-Sloan Atlas elliptical Petrosian analysis
 - S/N measurements
 - Binned spectra
 - Stellar kinematics and stellar-continuum models
 - Emission-line properties and models
 - Spectral-index measurements

It is important that you understand the procedures and algorithms used by the DAP to provide the data you're interested in if you are to trust their scientific usage.

'''Please consult the DAP [https://arxiv.org/pdf/1901.00856.pdf Overview] and [https://arxiv.org/pdf/1901.00866.pdf Emission-line Modeling] papers.'''

''Current release:'' '''MPL-8'''[[br]]
''DRP Version:'' '''v2_5_3'''[[br]]
''DAP Version:'' '''[https://github.com/sdss/mangadap/releases/tag/2.3.0 2.3.0]'''

{{{#!th
[[span(style=color: black, '''NOTE''')]]
}}}
{{{#!td
This page is meant to provide a useful guide for how to get started with the DAP output data.  ''If anything is unclear to you then this page is not doing its job! ''

Please '''[#QuestionsComments contact me]''' if you have any questions.

Also, if you have helpful code snippets that you'd like to share, please do!  You can put your code on the [https://trac.sdss.org/wiki/MANGA/code_share informal code wiki], or you can send the code directly to [#QuestionsComments Kyle] or the [https://github.com/sdss/marvin Marvin team] and we will include your code for general use in one of those repos.  You can also contribute to the [https://github.com/sdss/mangadap DAP] and/or [https://github.com/sdss/marvin Marvin] directly via [https://github.com/ GitHub].
}}}

## Directory structure ##

See description [wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#DirectoryStructure here].

## {{{[type]}}} selection ##

The subdirectories at the top level of the directory structure are named after the analysis {{{[type]}}}.  The {{{[type]}}} directories are:
 - {{{SPX-MILESHC-MILESHC/}}} - Analysis of each individual spaxel; spaxels must have a valid continuum fit for an emission-line model to be fit
 - {{{VOR10-MILESHC-MILESHC/}}} - Analysis of spectra binned to S/N~10 using the Voronoi binning algorithm (Cappellari & Copin 2003) 
 - {{{HYB10-MILESHC-MILESHC/}}} - Stellar-continuum analysis of spectra binned to S/N~10 for the stellar kinematics (same as VOR10 approach); however, the emission-line measurements are performed on the individual spaxels.  See a description of the hybrid binning scheme [wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#HYBbinningscheme here] and [wiki:MANGA/TRM/TRM_ActiveDev/dap#Executedanalysissteps here].

For more detail about the analyses performed, see the [wiki:MANGA/TRM/TRM_ActiveDev/dap#Executedanalysissteps analysis steps] in the DAP overview section.

The type of analysis that you should use will depend on your science application:

### {{{SPX-MILESHC-MILESHC/}}} ###

These are useful for most science applications that can push to very low S/N.  They're also useful for characterizing the performance of the measurements toward the low S/N limit of the data.

{{{#!th
[[span(style=color: red, '''WARNING''')]]
}}}
{{{#!td
 * Spectra with g-band S/N < 1 will not have a stellar-continuum model or Gaussian emission-line model.
 * '''Please consult the DAP [https://arxiv.org/pdf/1901.00856.pdf Overview] and [https://arxiv.org/pdf/1901.00866.pdf Emission-line Modeling] papers''' for usage guidelines and limitations of the data.
}}}

### {{{VOR10-MILESHC-MILESHC/}}} ###

These data are geared toward stellar kinematics, where we've Voronoi binned the spectra to a g-band S/N of approximately 10.

{{{#!th
[[span(style=color: red, '''WARNING''')]]
}}}
{{{#!td
 * No spectrum with a g-band S/N < 1 is included in any bin.
 * Voronoi binned spectra are just simple means of all the spectra in the bin.
 * The covariance in the datacube is propagated to the variance in the stacked spectra.
 * The spectral resolution in each binned spectra is propagated from the PREDISP cube provided by the DRP, similar to the formalism explained [[https://sdss-mangadap.readthedocs.io/en/latest/mangadap.drpfits.html#mangadap.drpfits.DRPFits.instrumental_dispersion_plane|here]].
 * (Binned) Spectra with g-band S/N < 1 will not have a stellar-continuum model or Gaussian emission-line model.
 * '''Please consult the DAP [https://arxiv.org/pdf/1901.00856.pdf Overview] and [https://arxiv.org/pdf/1901.00866.pdf Emission-line Modeling] papers''' for usage guidelines and limitations of the data.
 * Because the binning is done based on the ''continuum'' S/N, this limits the emission-line science that can be done at low continuum S/N.
}}}

### {{{HYB10-MILESHC-MILESHC/}}} ###

These are the default files that most users will want to use.  We first Voronoi binned the spectra to a g-band S/N of approximately 10 to measure the stellar kinematics.  Then these bins are deconstructed to fit the emission lines.  Coding/timescale issues prevented from providing spectral indices for this approach in this release.

{{{#!th
[[span(style=color: red, '''WARNING''')]]
}}}
{{{#!td
 * No spectrum with a g-band S/N < 1 is included in any bin.
 * Voronoi binned spectra are just simple means of all the spectra in the bin.
 * The covariance in the datacube is propagated to the variance in the stacked spectra.
 * The spectral resolution in each binned spectra is propagated from the PREDISP cube provided by the DRP, similar to the formalism explained [[https://sdss-mangadap.readthedocs.io/en/latest/mangadap.drpfits.html#mangadap.drpfits.DRPFits.instrumental_dispersion_plane|here]].
 * (Binned) Spectra with g-band S/N < 1 will not have a stellar-continuum model.
 * '''Please consult the DAP [https://arxiv.org/pdf/1901.00856.pdf Overview] and [https://arxiv.org/pdf/1901.00866.pdf Emission-line Modeling] papers''' for usage guidelines and limitations of the data.
 * All spectra with 80% valid pixels will have a combined emission-line+stellar-continuum model, where the stellar kinematics have been fixed by the fits to the binned spectra.
 * This is the only file where the BINIDs are different for the emission-line properties.
}}}


## Output files ##

The primary output files are located at:
 * [https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-8]/[type]/[plate]/[ifudesign].

There are two main output files for each observation (plate-ifudesign combination):
 * {{{manga-[PLATE]-[IFUDESIGN]-MAPS-[type].fits.gz}}}, see [#OutputMAPSfiles here]: 2D "maps" (i.e., images) of DAP measured properties
 * {{{manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[type].fits.gz}}}, see [#OutputLOGCUBEfiles here]: 3D data cubes with the binned and best-fitting-model spectra

The datacubes produced by the DAP have the same shape as the DRP datacube, and the DAP maps have the same on-sky pixel size as a single wavelength channel in the DRP datacubes.  This is meant to ease associating the DRP input and DAP output products.

Examples are given below for how to interact with the two main output files using python.  However, you are '''strongly encouraged''' to [http://sdss-marvin.readthedocs.io/en/stable/installation.html install Marvin] and use it to interact with the data.

### {{{HDUCLASS}}} ###

See [wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel/#HDUCLASS here]

### Output {{{MAPS}}} files ###

The {{{MAPS}}} files are the primary output file from the DAP; please see the detailed description of their format '''[wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#DAPMAPSfile here]'''.

In brief, the file contains 2D "maps" (i.e., images) of DAP measured properties.  Most properties are provided in groups of three fits extensions:
  1. {{{[property]}}}: the measurement value,
  2. {{{[property]_IVAR}}}: the measurement uncertainty stored as the inverse variance, and
  3. {{{[property]_MASK}}}: a corresponding bit mask for each spaxel.
The list of measured properties is described by the '''[wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#DAPMAPSfile data model]'''.

The headers of each extension provides the astrometric World Coordinate System (WCS) and should exactly match that of the DRP output {{{LOGCUBE}}} files (apart from the wavelength coordinate).

Many properties have multiple "species" or channels associated with them.  The identifying name of each mapped property is provided in the header; e.g., the emission-line channels are listed '''[wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel/#emlines here]'''.  In python, you can create a dictionary of items in each channel as follows:

{{{
#!python 

# Declare a function that creates a dictionary for the columns in the
# multi-channel extensions
def channel_dictionary(hdu, ext):
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict
}}}

or if you've installed the DAP, e.g.:

{{{
#!python 

from mangadap.util.fileio import channel_dictionary
from astropy.io import fits

hdu = fits.open('mangadap-7495-12704-MAPS-SPX-MILESHC-MILESHC.fits.gz')
emlc = channel_dictionary(hdu, 'EMLINE_GFLUX')
}}}

'''If possible''' with your software package of choice, please always select the extension and channel based on its '''name''', ''not'' its extension or channel number.  An example of how to do this using python is provided [#pythonusageexample below].  The ordering of, e.g., the emission lines in the relevant extensions has changed between different MPLs and may change again.

#### Necessary corrections ####

(''This information is a repeat of what is provided [wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#MAPcorrections here]'').

Note that the stellar and gas velocity dispersions '''must be corrected for instrumental resolution effects''' to obtain the astrophysical Doppler broadening.

The corrected gas velocity dispersion is:

 sigma_gas_corr = sqrt( {{{EMLINE_GSIGMA}}}^2^ - {{{EMLINE_INSTSIGMA}}}^2^ )
 
The corrected stellar velocity dispersion is:

 sigma_star_corr = sqrt( {{{STELLAR_SIGMA}}}^2^ - {{{STELLAR_SIGMACORR}}}^2^ )
 
'''In both cases''', beware of imaginary numbers.  That is, when the correction is larger than the provided value, the above equations result in taking the sqrt of a negative number.  Stellar velocity dispersions are provided for two approaches to the calculation

A nominal correction is calculated using the quadrature difference between the instrumental dispersion of the template and galaxy spectra over the fitted wavelenghth range.  This is the correction provided in MPL-5, MPL-7/DR15.  '''Please consult the DAP [https://arxiv.org/pdf/1901.00856.pdf Overview] paper.'''  In MPL-8, we also provide a correction based on a fit of the optimal template with and without the resolution matched to the MaNGA data; see more detail [wiki:MANGA/TRM/TRM_ActiveDev/dap#Executedanalysissteps here].

Also, velocity-dispersion corrections are provided for the spectral indices.  To apply the corrections, you have to know the unit of each index.  For angstrom units:

 specindex_ang_corr = {{{SPECINDEX}}} * {{{SPECINDEX_CORR}}}

and for magnitude units:

 specindex_mag_corr = {{{SPECINDEX}}} + {{{SPECINDEX_CORR}}}

An example of how to apply these corrections in detail is given below.

#### python usage example ####

With the {{{MAPS}}} FITS file, you should be able to extract DAP maps output using any fits reader.  ''Please let us know if you run into any problems! ''

For example, here is a code snippet that will plot the H-alpha flux map, stellar velocity field, the corrected stellar velocity dispersion field, and the corrected H-beta index map for '''[https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-8/HYB10-MILESHC-MILESHC/8138/12704/manga-8138-12704-MAPS-HYB10-MILESHC-MILESHC.fits.gz this MAPS]''' file using python:

{{{
#!python 

# Imports
from astropy.io import fits
from matplotlib import pyplot
import numpy

def channel_dictionary(hdu, ext):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict


def channel_units(hdu, ext):
    """
    Construct an array with the channel units.
    """
    nchannels = 1 if len(hdu[ext].data.shape) == 2 else hdu[ext].data.shape[0]
    channel_units = numpy.empty(nchannels, dtype=object)
    for k, v in hdu[ext].header.items():
        if k[0] == 'U':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_units[i] = v.strip()
    return channel_units


def apply_index_dispersion_correction(indx, indxcorr, unit):
    """
    Apply a set of dispersion corrections.
    """
    if unit not in [ 'ang', 'mag' ]:
        raise ValueError('Unit must be mag or ang.')
    return indx * indxcorr if unit == 'ang' else indx + indxcorr


# Open the fits file
hdu = fits.open('manga-8138-12704-MAPS-HYB10-MILESHC-MILESHC.fits.gz')

# Build a dictionary with the emission-line and spectral-index channel names to ease selection and get the spectral-index units
emlc = channel_dictionary(hdu, 'EMLINE_GFLUX')
spic = channel_dictionary(hdu, 'SPECINDEX')
spiu = channel_units(hdu, 'SPECINDEX')

# Show the Gaussian-fitted H-alpha flux map
mask_ext = hdu['EMLINE_GFLUX'].header['QUALDATA']
halpha_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['Ha-6564'],:,:], mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:] > 0)

pyplot.imshow(halpha_flux, origin='lower', interpolation='nearest', cmap='inferno')
pyplot.colorbar()
pyplot.show()

# Show the stellar velocity field
mask_ext = hdu['STELLAR_VEL'].header['QUALDATA']
stellar_vfield = numpy.ma.MaskedArray(hdu['STELLAR_VEL'].data, mask=hdu[mask_ext].data > 0)

pyplot.imshow(stellar_vfield, origin='lower', interpolation='nearest', vmin=-300, vmax=300, cmap='RdBu_r')
pyplot.colorbar()
pyplot.show()

# Show the corrected stellar velocity dispersion field
mask_ext = hdu['STELLAR_SIGMA'].header['QUALDATA']
stellar_sfield_sqr = numpy.ma.MaskedArray(numpy.square(hdu['STELLAR_SIGMA'].data) - numpy.square(hdu['STELLAR_SIGMACORR'].data[0,:,:]),
                                          mask=hdu[mask_ext].data > 0)
# WARNING: This will ignore any data where the correction is larger than the measurement
stellar_sfield = numpy.ma.sqrt(stellar_sfield_sqr)

pyplot.imshow(stellar_sfield, origin='lower', interpolation='nearest', cmap='viridis')
pyplot.colorbar()
pyplot.show()

# Show the corrected H-beta index map
mask_ext = hdu['SPECINDEX'].header['QUALDATA']
hbeta_index_raw = numpy.ma.MaskedArray(hdu['SPECINDEX'].data[spic['Hb'],:,:], mask=hdu[mask_ext].data[spic['Hb'],:,:] > 0)
hbeta_index = apply_index_dispersion_correction(hbeta_index_raw, hdu['SPECINDEX_CORR'].data[spic['Hb'],:,:], spiu[spic['Hb']])

pyplot.imshow(hbeta_index, origin='lower', interpolation='nearest', cmap='inferno')
pyplot.colorbar()
pyplot.show()

}}}

### Output {{{LOGCUBE}}} files ###

The {{{LOGCUBE}}} files provide the binned spectra and the best-fitting model spectrum for each spectrum that was successfully fit.  Please see the detailed description of their format '''[wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#DAPLOGCUBEfile here]'''.

These files are useful for detailed assessments of the model parameters because they allow you to return to the spectra and compare the model against the data.  As described by the DAP [https://arxiv.org/pdf/1901.00856.pdf Overview] paper, the DAP fits the spectra in two stages, one to get the stellar kinematics and the second to  determine the emission-line properties.  The emission-line module (used for all binning schemes) fits both the stellar continuum and the emission lines at the same time, where the stellar kinematics are fixed by the first fit.  The stellar-continuum models from the first fit are provided in the {{{STELLAR}}} extension; to get the stellar continuum determined during the emission-line modeling, you calculate:

stellar_continuum (from emission-line modeling module) = {{{MODEL}}} - {{{EMLINE}}}

{{{#!th
[[span(style=color: red, '''WARNING''')]]
}}}
{{{#!td
In the HYB binning case the binned spectra provided in the {{{LOGCUBE}}} files are from the Voronoi binning step.  However, the emission-line models are fit to the ''individual spaxels''.  So:
 * The stellar continuum fits from the first iteration, in the {{{STELLAR}}} extension, should be compared to the Voronoi binned spectra in the file, but
 * the best-fitting model spectra in the {{{MODEL}}} extension should be compared to the individual spectra from the DRP LOGCUBE file!
}}}

An example of how to plot the model cube data using python is provided below.

#### python usage example ####

With the {{{LOGCUBE}}} FITS file, you should be able to extract the binned spectra and best-fitting models produced by the DAP using any fits reading software.  ''Please let us know if you run into any problems! ''

For example, here is a code snippet that plots the highest S/N spectrum, the full model, the residuals, the model stellar continuum, and the model emission-line spectrum using '''[https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-8/VOR10-MILESHC-MILESHC/8138/12704/manga-8138-12704-MAPS-VOR10-MILESHC-MILESHC.fits.gz this MAPS]''' file and '''[https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-8/SPX-MILESHC-MILESHC/8138/12704/manga-8138-12704-LOGCUBE-VOR10-MILESHC-MILESHC.fits.gz this LOGCUBE]''' file using python:
 
{{{
#!python 

# Imports
from astropy.io import fits
from matplotlib import pyplot
import numpy

# This is a bitmask handling object from the DAP source code
from mangadap.dapfits import DAPCubeBitMask

# Open the fits file
hdu_maps = fits.open('manga-8138-12704-MAPS-SPX-MILESHC-MILESHC.fits.gz')
hdu_cube = fits.open('manga-8138-12704-LOGCUBE-SPX-MILESHC-MILESHC.fits.gz')

# Get the S/N per bin from the MAPS file
snr = numpy.ma.MaskedArray(hdu_maps['BIN_SNR'].data, mask=hdu_maps['BINID'].data[0,:,:] < 0)

# Select the bin/spaxel with the highest S/N
k = numpy.ma.argmax(snr.ravel())
n = hdu_maps['BIN_SNR'].data.shape[0] # Number of pixels in X and Y
# Get the pixel coordinate
j = k//n
i = k - j*n

# Declare the bitmask object to mask selected pixels
bm = DAPCubeBitMask()
wave = hdu_cube['WAVE'].data
flux = numpy.ma.MaskedArray(hdu_cube['FLUX'].data[:,j,i],
                            mask=bm.flagged(hdu_cube['MASK'].data[:,j,i],
				[ 'IGNORED', 'FLUXINVALID', 'IVARINVALID', 'ARTIFACT' ]))
model = numpy.ma.MaskedArray(hdu_cube['MODEL'].data[:,j,i],
                             mask=bm.flagged(hdu_cube['MODEL_MASK'].data[:,j,i], 'FITIGNORED'))
stellarcontinuum = numpy.ma.MaskedArray(
                        hdu_cube['MODEL'].data[:,j,i] - hdu_cube['EMLINE'].data[:,j,i],
                             mask=bm.flagged(hdu_cube['MODEL_MASK'].data[:,j,i], 'FITIGNORED'))
emlines = numpy.ma.MaskedArray(hdu_cube['EMLINE'].data[:,j,i],
                               mask=bm.flagged(hdu_cube['MODEL_MASK'].data[:,j,i],
                                               'ELIGNORED'))
resid = flux-model-0.5

pyplot.step(wave, flux, where='mid', color='k', lw=0.5)
pyplot.plot(wave, model, color='r', lw=1)
pyplot.plot(wave, stellarcontinuum, color='g', lw=1)
pyplot.plot(wave, emlines, color='b', lw=1)
pyplot.step(wave, resid, where='mid', color='0.5', lw=0.5)
pyplot.show()
}}}

## Using the pixel/spaxel masks ##

The maskbits for the DAP data are described '''[wiki:MANGA/TRM/TRM_MPL-8/DAPMetaData#Maskbits here]'''.  In particular, be aware of the DONOTUSE and the [wiki:MANGA/TRM/TRM_MPL-8/DAPMetaData#TheUNRELIABLEflag UNRELIABLE] flags for the {{{MAPS}}} files.

The 2d {{{MAPS}}} file pixel mask is [wiki:MANGA/TRM/TRM_MPL-8/DAPMetaData/#MANGA_DAPPIXMASK MANGA_DAPPIXMASK]

The 3d {{{LOGCUBE}}} file spaxel mask is [wiki:MANGA/TRM/TRM_MPL-8/DAPMetaData/#MANGA_DAPSPECMASK MANGA_DAPSPECMASK]

In all cases, the DAP has a convenience class that allows a user to quickly determine if any mask value is flagged with a certain value.  For example:
{{{
#!python 

# Imports
import os
from astropy.io import fits
from mangadap.util.bitmask import BitMask

# Define the path to the IDLUTILS maskbits file
sdssMaskbits = os.path.join(os.environ['IDLUTILS_DIR'], 'data', 'sdss', 'sdssMaskbits.par')

# Instantiate the BitMask object
bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_DAPQUAL')

# Read a DAP file
hdu = fits.open('manga-8138-12704-MAPS-SPX-MILESHC-MILESHC.fits.gz')

# Check if the file is critical and print the result
dap_file_is_critical = bm.flagged(hdu['PRIMARY'].header['DAPQUAL'], flag='CRITICAL')
print('This DAP file {0} flagged as CRITICAL.'.format('is' if dap_file_is_critical else 'is not'))
}}}

There are also a predefined set of derived [https://sdss-mangadap.readthedocs.io/en/latest/mangadap.util.bitmask.html#mangadap.util.bitmask.BitMask BitMask] classes that will use
the {{{$IDLUTILS_DIR}}} environmental variable to instantiate the
bitmask object.  For example:
{{{
#!python 

#Imports
import numpy
from astropy.io import fits
from mangadap.drpfits import DRPFitsBitMask

# Instantiate the BitMask object
#  - This requires an $IDLUTILS_DIR environmental variable
bm = DRPFitsBitMask()

# Read a DRP file
hdu = fits.open('manga-8138-12704-LOGCUBE.fits.gz')

# Find the number of pixels flagged as DONOTUSE or FORESTAR
indx = bm.flagged(hdu['MASK'].data, flag=['DONOTUSE', 'FORESTAR']) 
print('This DRP file has {0}/{1} pixels flagged as either DONOTUSE or FORESTAR.'.format(numpy.sum(indx), numpy.prod(indx.shape)))
}}}

For an example using IDL, see [wiki:MANGA/TRM/TRM_MPL-8/metadata#Bitmasks here].

'''There is also now a [http://sdss-marvin.readthedocs.io/en/stable/utils/maskbit.html Maskbits class] in Marvin!'''

## Using the BINID extension ##

The BINID extension has 5 channels.  They provide the IDs of spaxels associated with
 - 0. each binned spectrum.  Any spaxel with BINID=-1 as not included in any bin.
 - 1. any binned spectrum with an attempted stellar kinematics fit.
 - 2. any binned spectrum with emission-line moment measurements.
 - 3. any binned spectrum with an attempted emission-line fit.
 - 4. any binned spectrum with spectral-index measurements.

In any of these channels, you can obtain the unique bin numbers using {{{numpy.unique(bin_indx.ravel())[1:]}}} in python; the selection of all but the first array element is just provided all the numbers without the -1 for invalid spaxels.  If you're working with anything but the SPX binning, you'll want to extract the unique spectra and/or maps values.  You can do that by finding the indices of the unique bins, like this:
{{{
#!python 

unique_bins, unique_indices = tuple(map(lambda x : x[1:], numpy.unique(bin_indx.ravel(), return_index=True)))
}}}

Here's a worked example where I pull out the unique stellar velocities and produce a scatter plot of the x and y positions of the luminosity-weighted bin centers and color them by the measure stellar velocity.
{{{
#!python 

#Imports
import numpy
from astropy.io import fits
from matplotlib import pyplot

# Get the unique bins and indices
def unique_bins(bin_indx, return_index=False):
    """
    Get the unique bins andthe indices of the unique bins in the
    flattened spatial dimension, ignoring the bins with indices of -1.

    Same as mangadap.util.fitsutil.DAPFitsUtil.unique_bins.
    """
    return tuple(map(lambda x : x[1:], numpy.unique(bin_indx.ravel(), return_index=True))) \
                if return_index else numpy.unique(bin_indx.ravel())[1:]

# Read a DAP MAPS file
hdu = fits.open('manga-8138-12704-MAPS-HYB10-MILESHC-MILESHC.fits.gz')

# Get the unique indices of the stellar-kinematics bins
ubins, uindx = unique_bins(hdu['BINID'].data[1,:,:], return_index=True)

# Get the x and y coordinates and the stellar velocities
x = hdu['BIN_LWSKYCOO'].data[0,:,:].ravel()[uindx]
y = hdu['BIN_LWSKYCOO'].data[1,:,:].ravel()[uindx]
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
}}}

## Download via rsync ##

Here are some useful rsync commands that will

 * Grab all {{{MAPS}}}:

{{{
rsync -avz --include "*/" --include "*MAPS*fits.gz" --exclude "*" rsync://sdss@dtn01.sdss.utah.edu/sas/mangawork/manga/spectro/analysis/MPL-8/ path/to/mangadap/data/MPL-8/
}}}

 * Grab all {{{MAPS}}} for a single {{{[type]}}}:
  
{{{
rsync -avz --include "*/" --include "*MAPS*fits.gz" --exclude "*" rsync://sdss@dtn01.sdss.utah.edu/sas/mangawork/manga/spectro/analysis/MPL-8/[type] path/to/mangadap/data/MPL-8/[type]
}}}

 * Grab all {{{LOGCUBE}}}s:

{{{
rsync -avz --include "*/" --include "*LOGCUBE*fits.gz" --exclude "*" rsync://sdss@dtn01.sdss.utah.edu/sas/mangawork/manga/spectro/analysis/MPL-8/ path/to/mangadap/data/MPL-8/
}}}

 * Grab all {{{LOGCUBE}}}s for a single {{{[type]}}}:
  
{{{
rsync -avz --include "*/" --include "*LOGCUBE*fits.gz" --exclude "*" rsync://sdss@dtn01.sdss.utah.edu/sas/mangawork/manga/spectro/analysis/MPL-8/[type] path/to/mangadap/data/MPL-8/[type]
}}}

{{{#!th
[[span(style=color: black, '''NOTE''')]]
}}}
{{{#!td

 - You should change {{{path/to/mangadap/data/MPL-8/}}} to the directory where you want to hold the data on your local directory structure.
 - If you're using Marvin and you want to make sure you have the files locally, make sure that you make this path match where Marvin will look for the files!  See [http://sdss-marvin.readthedocs.io/en/stable/installation.html#local-sas-directory-structure here].
}}}

### Required disk space ###

For all {{{MAPS}}}:
||= {{{[type]}}} =||= '''Size (GB)''' =|| 
|| {{{SPX-MILESHC-MILESHC}}} || 15.1 ||
|| {{{VOR10-MILESHC-MILESHC}}} || 6.4 ||
|| {{{HYB10-MILESHC-MILESHC}}} || 16.1 ||

For all {{{LOGCUBE}}}s:
||= {{{[type]}}} =||= '''Size (GB)''' =||
|| {{{SPX-MILESHC-MILESHC}}} || 577.0 ||
|| {{{VOR10-MILESHC-MILESHC}}} || 201.0 ||
|| {{{HYB10-MILESHC-MILESHC}}} || 358.7 ||

## Product Certifications ##

For the purposes of development, we have adopted a rating system for each provided quantity from the DAP:

  - '''Unrated''': Little or no testing applied.  Use at your own risk.
  - '''B-rated certification''': Only minimal checks undertaken.  Products may be used for exploratory purposes, but users should test their science applications.
  - '''A-rated certification''': "Science ready" products have undergone basic verification that the algorithms have been optimized at the basic level and the limitations of the algorithms understood such that the data is suitably flagged.  Further tests may be needed for specific applications, especially to understand the limits of flagged data.
  - '''AAA-rated certification''': A thorough and comprehensive set of tests has been applied and the results meet or exceed survey requirements.

[[span(style=color: red, '''PENDING DPC APPROVAL''')]]

For MPL-8, we are recommending the following certifications by the data products committee (DPC):

'''Global DAP performance:''' A-rated[[br]]
A variety of tests indicate that catastrophic failures are limited to a few percent and are suitably identified.

'''Stellar Kinematics:''' A-rated[[br]]
The stellar kinematics have been significantly tested, and our recommended usage guidelines are provided [https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-8/DAPDataModel#Stellarvelocitydispersions here].  However, the {{{UNRELIABLE}}} flag has been turned off, and the flagging is even more limited than the MPL-5 data.  This is because of worries over biasing averaged kinematics when one simply removes the low dispersion, low S/N measurements.  Further discussion is provided in the DAP paper (for now see [https://trac.sdss.org/attachment/wiki/MANGA/MaNGA/Meetings/Mexico_Dec2017/westfall_stellarkin.pdf this] presentation from the Mexico team meeting).

'''Emission Line Measurements:'''
 Strong-line fluxes and EWs: A-rated[[br]]
 Weak-line fluxes and EWs: B-rated[[br]]
 Kinematics: B-rated[[br]]

The emission-line module is still relatively new, and it is substantially different from the one in MPL-5.  We have checked the fluxes of a few of the strong lines (H-alpha, H-beta, OIII), but much more testing is needed, particularly for the weak lines.  The velocities of all lines are tied, which may lead to biases in the fits given the different ionization potentials and geometries expected between the Balmer and higher ionization lines.  These biases could be present in effectively all the reported emission-line parameters.  The velocity dispersions of most lines are independent, and the line widths could have significant systematic errors toward low S/N.  Further testing and characterization of the line-fitting limitations will continue toward the writing of the DAP paper.

'''Spectral indices:'''
 D4000 and Dn4000: A-rated[[br]]
 Other indices: B-rated

The algorithm used for the spectral-index measurements is largely unchanged since MPL-5.  The measurements have been compared to those provided in the DR14 Firefly VAC (see [http://www.sdss.org/dr14/manga/manga-data/manga-firefly-value-added-catalog/ here]), and they are fully consistent with one another with a couple of caveats:
  - The index measurements are provided as ''uncorrected'' by the DAP.  I.e., the DAP index measurements do not account for the effect of the velocity dispersion on the index.  Both the Firefly VAC and DAP provide the "velocity-dispersion corrections"; when using the DAP, these corrections '''have to be applied to the data by the user''', as described [wiki:MANGA/TRM/TRM_MPL-8/DAPDataModel#MAPcorrections here].
  - The errors provided by the DAP have been shown to be smaller than the errors reported by the Firefly VAC, and the latter should be more formally correct.  The systematic overestimation of the errors is true of ''all'' the indices, even though the D4000 and Dn4000 indices have been A-rated.

## Questions?  Comments? ##

Please contact [mailto:westfall@ucolick.org Kyle Westfall] if you have any problems, comments, suggestions, bug reports, etc.

