import os
import time
import argparse
import numpy

from scipy.special import erf
from scipy import interpolate

from matplotlib import pyplot, ticker, rc, colors, cm, colorbar, image

from astropy.io import fits

from mangadap.config import defaults
from mangadap.proc.util import growth_lim
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel

from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.util.filter import BoxcarFilter
from mangadap.util.sampling import spectral_coordinate_step
from mangadap.util.mapping import map_extent, map_beam_patch
from mangadap.util.fileio import channel_dictionary

from mangadap.scripts import scriptbase


def init_ax(fig, pos, facecolor=None, grid=False):
    ax = fig.add_axes(pos, facecolor=facecolor)
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    if grid:
        ax.grid(True, which='major', color='0.8', zorder=1, linestyle='-', lw=0.5)
    return ax


def init_image_ax(fig, pos):
    ax = fig.add_axes(pos, facecolor='0.9')
    ax.minorticks_on()
    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6, direction='in')
    ax.tick_params(which='minor', length=3, direction='in')
    return ax


def masked_imshow(fig, ax, cax, data, extent=None, norm=None, vmin=None, vmax=None, cmap='viridis',
                  zorder=0, cbformat=None, subs=None):
    if numpy.sum(data.mask) != numpy.prod(data.shape):
        if norm is None:
            img = ax.imshow(data, origin='lower', interpolation='nearest', extent=extent,
                            vmin=vmin, vmax=vmax, cmap=cmap, zorder=zorder)
        else:
            img = ax.imshow(data, origin='lower', interpolation='nearest', extent=extent,
                            norm=norm, cmap=cmap, zorder=zorder)

        cb = fig.colorbar(img, cax=cax) if cbformat is None else \
                    fig.colorbar(img, cax=cax, format=cbformat)
        if subs is not None:
            cb.locator = ticker.LogLocator(base=10, subs=(1.,2.,4.,))
            cb.update_ticks()
    else:
        _norm = colors.Normalize(vmin=vmin, vmax=vmax) if norm is None else norm
        cb = colorbar.ColorbarBase(cax, cmap=cm.get_cmap(cmap), norm=_norm)


def gmr_data(drp_file):
    if not os.path.isfile(drp_file):
        raise FileNotFoundError('{0} does not exist!'.format(drp_file))

    hdu = fits.open(drp_file)
    gmr_map = -2.5*numpy.ma.log10(numpy.ma.MaskedArray(hdu['GIMG'].data,
                                                       mask=numpy.invert(hdu['GIMG'].data>0))
                                    / numpy.ma.MaskedArray(hdu['RIMG'].data,
                                                           mask=numpy.invert(hdu['RIMG'].data>0)))
    hdu.close()
    del hdu
    return gmr_map


def get_stellar_continuum_data(dapmaps, dapcube):

    if not os.path.isfile(dapmaps):
        raise FileNotFoundError('{0} does not exist.'.format(dapmaps))
    if not os.path.isfile(dapcube):
        raise FileNotFoundError('{0} does not exist.'.format(dapcube))

    # Get the data
    cube_hdu = fits.open(dapcube)

    # Find the unique bins and get rid of -1 bins
    binid = cube_hdu['BINID'].data[1,:,:].ravel()
#    uniq, indx = map(lambda x: x[1:], numpy.unique(binid, return_index=True))
    uniq, indx = numpy.unique(binid, return_index=True)
    if uniq[0] == -1:
        uniq = uniq[1:]
        indx = indx[1:]

    # Get the unique fluxes, errors, and models
    npix = cube_hdu['FLUX'].data.shape[0]
    flux = numpy.ma.MaskedArray(cube_hdu['FLUX'].data.reshape(npix,-1),
                                mask=cube_hdu['MASK'].data.reshape(npix,-1) > 0)[:,indx]
    error = numpy.ma.power(numpy.ma.MaskedArray(cube_hdu['IVAR'].data.reshape(npix,-1),
                                                mask=cube_hdu['MASK'].data.reshape(npix,-1) > 0),
                           -0.5)[:,indx]
    model = numpy.ma.MaskedArray(cube_hdu['STELLAR'].data.reshape(npix,-1),
                                 mask=cube_hdu['STELLAR_MASK'].data.reshape(npix,-1) > 0)[:,indx]

    # Get the associated S/N and fit metrics
    maps_hdu = fits.open(dapmaps)
    snrg = maps_hdu['BIN_SNR'].data.ravel()[indx]
    rms = maps_hdu['STELLAR_FOM'].data[0,:,:].ravel()[indx]
    frms = maps_hdu['STELLAR_FOM'].data[1,:,:].ravel()[indx]
    rchi2 = maps_hdu['STELLAR_FOM'].data[2,:,:].ravel()[indx]

    return cube_hdu['WAVE'].data, flux, error, model, snrg, rms, frms, rchi2


def get_stellar_continuum_fom_maps(drpcube, dapmaps):

    if not os.path.isfile(dapmaps):
        raise FileNotFoundError('{0} does not exist.'.format(dapmaps))

    # Get the associated S/N and fit metrics
    maps_hdu = fits.open(dapmaps)
    snrg = numpy.ma.MaskedArray(maps_hdu['BIN_SNR'].data,
                                mask=numpy.invert(maps_hdu['BIN_SNR'].data > 0))
    mask = maps_hdu['STELLAR_VEL_MASK'].data > 0
    rms = numpy.ma.MaskedArray(maps_hdu['STELLAR_FOM'].data[0,:,:], mask=mask)
    frms = numpy.ma.MaskedArray(maps_hdu['STELLAR_FOM'].data[1,:,:], mask=mask)
    rchi2 = numpy.ma.MaskedArray(maps_hdu['STELLAR_FOM'].data[2,:,:], mask=mask)
    chi_grw = numpy.ma.MaskedArray(maps_hdu['STELLAR_FOM'].data[6:,:,:])
    chi_grw[0,mask] = numpy.ma.masked
    chi_grw[1,mask] = numpy.ma.masked
    chi_grw[2,mask] = numpy.ma.masked

    extent = map_extent(maps_hdu, 'SPX_MFLUX')

    return extent, snrg, gmr_data(drpcube), rms, frms, rchi2, chi_grw


def get_emission_line_data(drpcube, daptype, dapmaps, dapcube):
    if not os.path.isfile(dapmaps):
        raise FileNotFoundError('{0} does not exist.'.format(dapmaps))
    if not os.path.isfile(dapcube):
        raise FileNotFoundError('{0} does not exist.'.format(dapcube))

    # Get the data
    cube_hdu = fits.open(dapcube)
    maps_hdu = fits.open(dapmaps)

    # Get the observed spectra and errors
    if 'HYB' in daptype:
        if not os.path.isfile(drpcube):
            raise FileNotFoundError('{0} does not exist.'.format(drpcube))
        # Fits are to unbinned spectra from DRP datacube
        drp_hdu = fits.open(drpcube)
        npix = drp_hdu['FLUX'].data.shape[0]
        flux = numpy.ma.MaskedArray(drp_hdu['FLUX'].data.reshape(npix,-1),
                                    mask=drp_hdu['MASK'].data.reshape(npix,-1) > 0)
        error = numpy.ma.power(numpy.ma.MaskedArray(drp_hdu['IVAR'].data.reshape(npix,-1),
                                            mask=drp_hdu['MASK'].data.reshape(npix,-1) > 0), -0.5)
        drp_hdu.close()
        snrg = maps_hdu['SPX_SNR'].data.ravel()

    else:
        # Fits are to binned spectra from DAP datacube
        npix = cube_hdu['FLUX'].data.shape[0]
        flux = numpy.ma.MaskedArray(cube_hdu['FLUX'].data.reshape(npix,-1),
                                    mask=cube_hdu['MASK'].data.reshape(npix,-1) > 0)
        error = numpy.ma.power(numpy.ma.MaskedArray(cube_hdu['IVAR'].data.reshape(npix,-1),
                                                mask=cube_hdu['MASK'].data.reshape(npix,-1) > 0),
                               -0.5)
        snrg = maps_hdu['BIN_SNR'].data.ravel()
   
    # Find the unique bins and get rid of -1 bins
    binid = cube_hdu['BINID'].data[3,:,:].ravel()
#    uniq, indx = map(lambda x: x[1:], numpy.unique(binid, return_index=True))
    uniq, indx = numpy.unique(binid, return_index=True)
    if uniq[0] == -1:
        uniq = uniq[1:]
        indx = indx[1:]

    flux = flux[:,indx]
    error = error[:,indx]
    model = numpy.ma.MaskedArray(cube_hdu['MODEL'].data.reshape(npix,-1),
                                 mask=cube_hdu['MODEL_MASK'].data.reshape(npix,-1) > 0)[:,indx]

    # Get the associated S/N and fit metrics
    snrg = snrg[indx]
    rms = maps_hdu['EMLINE_FOM'].data[0,:,:].ravel()[indx]
    frms = maps_hdu['EMLINE_FOM'].data[1,:,:].ravel()[indx]
    rchi2 = maps_hdu['EMLINE_FOM'].data[2,:,:].ravel()[indx]

    return cube_hdu['WAVE'].data, flux, error, model, snrg, rms, frms, rchi2


def get_emission_line_fom_maps(dapmaps):

    maps_hdu = fits.open(dapmaps)
    eml = channel_dictionary(maps_hdu, 'EMLINE_GFLUX')
    anr = numpy.ma.MaskedArray(maps_hdu['EMLINE_GANR'].data[eml['Ha-6564'],...],
                               mask=maps_hdu['EMLINE_GFLUX_MASK'].data[eml['Ha-6564'],...] > 0)
    ha = numpy.ma.MaskedArray(maps_hdu['EMLINE_GFLUX'].data[eml['Ha-6564'],...],
                              mask=maps_hdu['EMLINE_GFLUX_MASK'].data[eml['Ha-6564'],...] > 0)
    hb = numpy.ma.MaskedArray(maps_hdu['EMLINE_GFLUX'].data[eml['Hb-4862'],...],
                              mask=maps_hdu['EMLINE_GFLUX_MASK'].data[eml['Hb-4862'],...] > 0)
    aob = numpy.ma.divide(ha, hb)

    # Get the associated S/N and fit metrics
    rms = numpy.ma.MaskedArray(maps_hdu['EMLINE_FOM'].data[0,:,:])
    rms[numpy.invert(rms > 0)] = numpy.ma.masked
    frms = numpy.ma.MaskedArray(maps_hdu['EMLINE_FOM'].data[1,:,:])
    frms[numpy.invert(frms > 0)] = numpy.ma.masked
    rchi2 = numpy.ma.MaskedArray(maps_hdu['EMLINE_FOM'].data[2,:,:])
    rchi2[numpy.invert(rchi2 > 0)] = numpy.ma.masked
    chi_grw = numpy.ma.MaskedArray(maps_hdu['EMLINE_FOM'].data[6:,:,:])
    chi_grw[numpy.invert(chi_grw > 0)] = numpy.ma.masked

    extent = map_extent(maps_hdu, 'SPX_MFLUX')

    return extent, anr, aob, rms, frms, rchi2, chi_grw


def fom_lambda(plt, ifu, wave, flux, error, model, fit='sc', wave_limits=None, ofile=None):

    if wave_limits is not None:
        indx = (wave < wave_limits[0]) | (wave > wave_limits[1])
        flux[indx,:] = numpy.ma.masked
        error[indx,:] = numpy.ma.masked
        model[indx,:] = numpy.ma.masked
    else:
        indx = numpy.arange(len(wave))[numpy.any(numpy.invert(model.mask), axis=1)]
        wave_limits = [wave[indx[0]], wave[indx[-1]]]

    nspec = flux.shape[1]
    resid = numpy.ma.absolute(flux-model)
    fresid = numpy.ma.absolute(numpy.ma.divide(resid, model))
    chi = numpy.ma.absolute(numpy.ma.divide(resid,error))
    snr = numpy.ma.divide(flux, error)

    mean_flux = numpy.ma.mean(flux, axis=1)
    mean_snr = numpy.ma.mean(snr, axis=1)
    mean_resid = numpy.ma.mean(resid, axis=1)
    mean_fresid = numpy.ma.mean(fresid, axis=1)
    mean_chi = numpy.ma.mean(chi, axis=1)

    bc = BoxcarFilter(100)
    resid_sm = bc.smooth(resid)
    fresid_sm = bc.smooth(fresid)
    chi_sm = bc.smooth(chi)

    spec_lim = [-10, nspec+10]

    flux_lim = numpy.exp(growth_lim(numpy.ma.log(mean_flux).compressed(), 0.99, fac=2.0))
    snr_lim = numpy.exp(growth_lim(numpy.ma.log(mean_snr).compressed(), 0.99, fac=2.0))
    resid_lim = numpy.exp(growth_lim(numpy.ma.log(mean_resid).compressed(), 0.99, fac=2.0))
    fresid_lim = numpy.exp(growth_lim(numpy.ma.log(mean_fresid).compressed(), 0.99, fac=2.0))
    chi_lim = numpy.exp(growth_lim(numpy.ma.log(mean_chi).compressed(), 0.99, fac=2.0, midpoint=0))

    l0 = numpy.log10(wave[0])
    dl = spectral_coordinate_step(wave, log=True)
    wave_bins = numpy.power(10, l0-dl/2 + dl*numpy.arange(len(wave)+1))
    spec_bins = 0.5+numpy.arange(nspec+1)

    font = { 'size' : 10 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.2, 0.87, 0.60, 0.12], facecolor='0.95', grid=True)
    ax.set_xlim(wave_limits)
    ax.set_ylim(flux_lim)
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.plot(wave, numpy.ma.mean(flux, axis=1), color='k', lw=0.5)
    ax.text(-0.15, 0.5, r'Flux', horizontalalignment='center',
            verticalalignment='center', rotation='vertical', transform=ax.transAxes)
    axt = ax.twinx()
    axt.set_xlim(wave_limits)
    axt.set_ylim(snr_lim)
    axt.set_yscale('log')
    axt.tick_params(axis='y', colors='0.5')
    axt.xaxis.set_major_formatter(ticker.NullFormatter())
    axt.plot(wave, numpy.ma.mean(snr, axis=1), color='0.5', lw=0.5)
    axt.text(1.1, 0.5, r'S/N', horizontalalignment='center',
            verticalalignment='center', rotation='vertical', transform=axt.transAxes,
            color='0.5')

    ax = init_ax(fig, [0.2, 0.75, 0.60, 0.12], facecolor='0.95', grid=True)
    ax.set_xlim(wave_limits)
    ax.set_ylim(resid_lim)
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
#    ax.plot(wave, numpy.ma.mean(resid_sm, axis=0), color='C3')
    ax.plot(wave, numpy.ma.mean(resid, axis=1), color='k', lw=0.5)
    ax.text(-0.15, 0.5, r'$|\Delta|$', horizontalalignment='center',
            verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    ax = init_ax(fig, [0.2, 0.63, 0.60, 0.12], facecolor='0.95', grid=True)
    ax.set_xlim(wave_limits)
    ax.set_ylim(fresid_lim)
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.plot(wave, numpy.ma.mean(fresid, axis=1), color='k', lw=0.5)
    ax.text(-0.15, 0.5, r'$|\Delta|/m$', horizontalalignment='center',
            verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    ax = init_ax(fig, [0.2, 0.51, 0.60, 0.12], facecolor='0.95', grid=True)
    ax.set_xlim(wave_limits)
    ax.set_ylim(chi_lim)
    ax.set_yscale('log')
    ax.plot(wave, numpy.ma.mean(chi, axis=1), color='k', lw=0.5)
    ax.text(0.5, -0.41, r'Wavelength ($\AA$)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.15, 0.5, r'$|\Delta|/\epsilon$', horizontalalignment='center',
            verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    label = '{0}-{1}'.format(plt, ifu)
    label += (' STELLAR' if fit == 'sc' else ' EMLINE')
    ax.text(1.23, -0.3, label, horizontalalignment='center',
            verticalalignment='bottom', rotation='vertical', transform=ax.transAxes,
            fontsize=22)

    ax = init_ax(fig, [0.07, 0.23, 0.30, 0.16], facecolor='0.95')
    ax.set_xlim(wave_limits)
    ax.set_ylim(spec_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    im = ax.pcolormesh(wave_bins, spec_bins, resid.T,
                       norm=colors.LogNorm(vmin=resid_lim[0], vmax=resid_lim[1]),
                       cmap='viridis', zorder=4, lw=0, rasterized=True)
    ax = init_ax(fig, [0.07, 0.07, 0.30, 0.16], facecolor='0.95')
    ax.set_xlim(wave_limits)
    ax.set_ylim(spec_lim)
    ax.pcolormesh(wave_bins, spec_bins, resid_sm.T,
                  norm=colors.LogNorm(vmin=resid_lim[0], vmax=resid_lim[1]),
                  cmap='viridis', zorder=4, lw=0, rasterized=True)
    cax = fig.add_axes([0.07, 0.395, 0.20, 0.01])
    cb = pyplot.colorbar(im, cax=cax, orientation='horizontal') #, format=FormatStrFormatter('%d'))
    cb.ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False,
                      labeltop=True)
    cax.text(1.05, 0.7, r'$|\Delta|$', horizontalalignment='left', verticalalignment='center',
            transform=cax.transAxes, fontsize=10)
    
    ax = init_ax(fig, [0.37, 0.23, 0.30, 0.16], facecolor='0.95')
    ax.set_xlim(wave_limits)
    ax.set_ylim(spec_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    im = ax.pcolormesh(wave_bins, spec_bins, fresid.T,
                       norm=colors.LogNorm(vmin=fresid_lim[0], vmax=fresid_lim[1]),
                       cmap='viridis', zorder=4, lw=0, rasterized=True)
    ax = init_ax(fig, [0.37, 0.07, 0.30, 0.16], facecolor='0.95')
    ax.set_xlim(wave_limits)
    ax.set_ylim(spec_lim)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.pcolormesh(wave_bins, spec_bins, fresid_sm.T,
                  norm=colors.LogNorm(vmin=fresid_lim[0], vmax=fresid_lim[1]),
                  cmap='viridis', zorder=4, lw=0, rasterized=True)
    ax.text(0.5, -0.30, r'Wavelength ($\AA$)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    cax = fig.add_axes([0.37, 0.395, 0.20, 0.01])
    cb = pyplot.colorbar(im, cax=cax, orientation='horizontal') #, format=FormatStrFormatter('%d'))
    cb.ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False,
                      labeltop=True)
    cax.text(1.05, 0.7, r'$|\Delta|/m$', horizontalalignment='left', verticalalignment='center',
            transform=cax.transAxes, fontsize=10)

    ax = init_ax(fig, [0.67, 0.23, 0.30, 0.16], facecolor='0.95')
    ax.set_xlim(wave_limits)
    ax.set_ylim(spec_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    im = ax.pcolormesh(wave_bins, spec_bins, chi.T,
                       norm=colors.LogNorm(vmin=chi_lim[0], vmax=chi_lim[1]),
                       cmap='RdBu_r', zorder=4, lw=0, rasterized=True)
    ax = init_ax(fig, [0.67, 0.07, 0.30, 0.16], facecolor='0.95')
    ax.set_xlim(wave_limits)
    ax.set_ylim(spec_lim)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.pcolormesh(wave_bins, spec_bins, chi_sm.T,
                  norm=colors.LogNorm(vmin=chi_lim[0], vmax=chi_lim[1]),
                  cmap='RdBu_r', zorder=4, lw=0, rasterized=True)
    cax = fig.add_axes([0.67, 0.395, 0.20, 0.01])
    cb = pyplot.colorbar(im, cax=cax, orientation='horizontal') #, format=FormatStrFormatter('%d'))
    cb.ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False,
                      labeltop=True)
    cax.text(1.05, 0.7, r'$|\Delta|/\epsilon$', horizontalalignment='left',
             verticalalignment='center', transform=cax.transAxes, fontsize=10)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def fom_growth(plt, ifu, flux, error, model, snrg, rms, frms, rchi2, fit='sc', ofile=None):

    nspec = flux.shape[1]
    resid = numpy.ma.absolute(flux-model)
    fresid = numpy.ma.divide(resid, model)
    chi = numpy.ma.divide(resid, error)

    snr_lim = numpy.exp(growth_lim(numpy.ma.log(snrg).compressed(), 0.95, fac=1.3))

    chi_lim = numpy.exp(growth_lim(numpy.ma.log(chi).compressed(), 0.8, fac=1.1))
    chi_lim[1] = numpy.ma.amax(chi)*2
    resid_lim = numpy.exp(growth_lim(numpy.ma.log(resid).compressed(), 0.9, fac=1.1))
    resid_lim[1] = numpy.ma.amax(resid)*2
    fresid_lim = numpy.exp(growth_lim(numpy.ma.log(fresid).compressed(), 0.9, fac=1.1))
    fresid_lim[1] = numpy.ma.amax(fresid)*2

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    rax = init_ax(fig, [0.07, 0.50, 0.30, 0.35], facecolor='0.9', grid=True)
    rax.set_xlim(resid_lim)
    rax.set_ylim([1e-4,2])
    rax.set_yscale('log')
    rax.set_xscale('log')
    rax.text(0.5, -0.17, r'$|\Delta|$', ha='center', va='center', transform=rax.transAxes)
    rax.text(0.1, 0.5, r'1-Growth', horizontalalignment='center',
             verticalalignment='center', transform=rax.transAxes, rotation='vertical')

    fax = init_ax(fig, [0.37, 0.50, 0.30, 0.35], facecolor='0.9', grid=True)
    fax.set_xlim(fresid_lim)
    fax.set_ylim([1e-4,2])
    fax.set_yscale('log')
    fax.set_xscale('log')
    fax.yaxis.set_major_formatter(ticker.NullFormatter())
    label = '{0}-{1}'.format(plt, ifu)
    label += (' STELLAR' if fit == 'sc' else ' EMLINE')
    fax.text(0.5, 1.3, label, horizontalalignment='center',
             verticalalignment='center', transform=fax.transAxes, fontsize=22)
    fax.text(0.5, -0.17, r'$|\Delta|/m$', ha='center', va='center', transform=fax.transAxes)

    cax = init_ax(fig, [0.67, 0.50, 0.30, 0.35], facecolor='0.9', grid=True)
    cax.set_xlim(chi_lim)
    cax.set_ylim([1e-4,2])
    cax.set_yscale('log')
    cax.set_xscale('log')
    cax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax.text(0.5, -0.17, r'$|\Delta|/\epsilon$', ha='center', va='center',
             transform=cax.transAxes)


    cmap = cm.get_cmap('viridis')
    cnorm = colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1])
    cbax = fig.add_axes([0.69, 0.86, 0.2, 0.01])

    cb = colorbar.ColorbarBase(cbax, norm=cnorm, cmap=cmap, orientation='horizontal')
    cb.ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False,
                      labeltop=True)
    cbax.text(1.05, 0.5, r'S/N$_g$', ha='left', va='center', transform=cbax.transAxes)

    snr_srt = numpy.argsort(snrg)

    for i in range(nspec):
        j = snr_srt[i]
        if numpy.all(chi.mask[:,j]):
            continue

        t = resid[:,j].compressed()
        s = numpy.argsort(t)
        rax.step(t[s], 1-numpy.arange(len(t)).astype(float)/len(t), where='pre',
                 color=cmap(cnorm(snrg[j])), lw=0.7)

        t = fresid[:,j].compressed()
        s = numpy.argsort(t)
        fax.step(t[s], 1-numpy.arange(len(t)).astype(float)/len(t), where='pre',
                 color=cmap(cnorm(snrg[j])), lw=0.7)

        t = chi[:,j].compressed()
        s = numpy.argsort(t)
        cax.step(t[s], 1-numpy.arange(len(t)).astype(float)/len(t), where='pre',
                 color=cmap(cnorm(snrg[j])), lw=0.7)

    x = numpy.linspace(0,5,100)
    g = (erf(x/numpy.sqrt(2)) - erf(-x/numpy.sqrt(2)))/2.
    interp = interpolate.interp1d(g, x)
    cax.plot(x,1-g, color='k')

    chi2_lim = numpy.exp(growth_lim(numpy.ma.log(rchi2).compressed(), 1.0, fac=2.0))
    rms_lim = numpy.exp(growth_lim(numpy.ma.log(rms).compressed(), 1.0, fac=2.0))
    frms_lim = numpy.exp(growth_lim(numpy.ma.log(frms).compressed(), 1.0, fac=2.0))
    snr_lim = numpy.exp(growth_lim(numpy.ma.log(snrg).compressed(), 1.0, fac=2.0))

    rax = init_ax(fig, [0.07, 0.15, 0.30, 0.25], facecolor='0.9', grid=True)
    rax.set_xlim(rms_lim)
    rax.set_ylim(snr_lim)
    rax.set_yscale('log')
    rax.set_xscale('log')
    rax.scatter(rms, snrg,
                marker='.', s=20, color='k', lw=0, alpha=0.5)
    rax.text(0.5, -0.25, r'$RMS$', ha='center', va='center', transform=rax.transAxes)
    rax.text(0.1, 0.5, r'S/N$_g$', horizontalalignment='center',
             verticalalignment='center', transform=rax.transAxes, rotation='vertical')

    fax = init_ax(fig, [0.37, 0.15, 0.30, 0.25], facecolor='0.9', grid=True)
    fax.set_xlim(frms_lim)
    fax.set_ylim(snr_lim)
    fax.set_yscale('log')
    fax.set_xscale('log')
    fax.yaxis.set_major_formatter(ticker.NullFormatter())
    fax.scatter(frms, snrg,
                marker='.', s=20, color='k', lw=0, alpha=0.5)
    fax.text(0.5, -0.25, r'$fRMS$', ha='center', va='center', transform=fax.transAxes)

    cax = init_ax(fig, [0.67, 0.15, 0.30, 0.25], facecolor='0.9', grid=True)
    cax.set_xlim(chi2_lim)
    cax.set_ylim(snr_lim)
    cax.set_yscale('log')
    cax.set_xscale('log')
    cax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax.scatter(rchi2, snrg,
                marker='.', s=20, color='k', lw=0, alpha=0.5)
    cax.text(0.5, -0.25, r'$\chi^2_\nu$', ha='center', va='center', transform=cax.transAxes)
    
    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def fom_maps(plt, ifu, image_file, snr, color, rms, frms, rchi2, chi_grw, extent=None,
             fit='sc', ofile=None):

    color_lim = numpy.exp(growth_lim(numpy.ma.log(color), 0.95, fac=1.1)) \
                    if fit != 'sc' else growth_lim(color, 0.95, fac=1.10)

    snr_lim = numpy.power(10., growth_lim(numpy.ma.log10(snr), 0.90, fac=1.05))
    rms_lim = numpy.power(10., growth_lim(numpy.ma.log10(rms), 0.90, fac=1.05))
    frm_lim = numpy.power(10., growth_lim(numpy.ma.log10(frms), 0.90, fac=1.05))

    chi_lim = [1/4, 4]
    c68_lim = [1/4, 4]
    c99_lim = [2.6/4, 2.6*4]
    crt_lim = [0.8, 2]

    left = 0.035
    bott = 0.08
    imwd = 0.24
    hbuf = 0.08
    cbuf = 0.005
    cbwd = 0.01
    vbuf = 0.03

    font = { 'size' : 6 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    i,j=0,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    if os.path.isfile(image_file):
        img = image.imread(image_file)
        ax.imshow(img)
        ax.text(0.95, 0.95, r'SDSS', color='white', ha='right', va='center', transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, 'No Image', horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=20)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    i,j=1,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, snr, extent=extent,
                  norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]), cmap='viridis',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'S/N$_g$' if fit == 'sc' else r'H$\alpha$ A/N', ha='right', va='center',
            transform=ax.transAxes, zorder=4)

    label = '{0}-{1}'.format(plt, ifu)
    label += (' STELLAR' if fit == 'sc' else ' EMLINE')
    ax.text(0.5, 1.1, label, horizontalalignment='center',
            verticalalignment='bottom', transform=ax.transAxes, fontsize=26)

    i,j=2,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    kwargs = {}
    if fit == 'sc':
        # Plotting the g-r color
        kwargs['vmin'] = color_lim[0]
        kwargs['vmax'] = color_lim[1]
    else:
        # Plotting H-alpha/H-beta flux ratio
        kwargs['norm'] = colors.LogNorm(vmin=color_lim[0], vmax=color_lim[1])
    masked_imshow(fig, ax, cax, color, extent=extent, cmap='RdBu_r', zorder=3, **kwargs)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'$g-r$' if fit == 'sc' else r'H$\alpha$/H$\beta$', ha='right',
            va='center', transform=ax.transAxes, zorder=4)

    i,j=0,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, rms, extent=extent,
                  norm=colors.LogNorm(vmin=rms_lim[0], vmax=rms_lim[1]), cmap='viridis',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'RMS', ha='right', va='center', transform=ax.transAxes, zorder=4)

    i,j=1,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, frms, extent=extent,
                  norm=colors.LogNorm(vmin=frm_lim[0], vmax=frm_lim[1]), cmap='viridis',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'fRMS', ha='right', va='center', transform=ax.transAxes, zorder=4)

    i,j=2,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, rchi2, extent=extent,
                  norm=colors.LogNorm(vmin=chi_lim[0], vmax=chi_lim[1]), cmap='RdBu_r',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'$\chi^2$', ha='right', va='center', transform=ax.transAxes, zorder=4)

    i,j=0,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, chi_grw[0], extent=extent,
                  norm=colors.LogNorm(vmin=c68_lim[0], vmax=c68_lim[1]), cmap='viridis',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'68% Growth $|\Delta|/\epsilon$', ha='right', va='center',
            transform=ax.transAxes, zorder=4)

    i,j=1,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, chi_grw[1], extent=extent,
                  norm=colors.LogNorm(vmin=c99_lim[0], vmax=c99_lim[1]), cmap='viridis',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'99% Growth $|\Delta|/\epsilon$', ha='right', va='center',
            transform=ax.transAxes, zorder=4)

    x = numpy.linspace(0,5,100)
    g = (erf(x/numpy.sqrt(2)) - erf(-x/numpy.sqrt(2)))/2.
    interp = interpolate.interp1d(g, x)
    m = interp([0.68, 0.99])
    i,j=2,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    masked_imshow(fig, ax, cax, numpy.ma.divide(chi_grw[1],chi_grw[0])*m[0]/m[1],
                  extent=extent, vmin=crt_lim[0], vmax=crt_lim[1], cmap='RdBu_r',
                  zorder=3)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.95, 0.95, r'99%/68% wrt Gaussian', ha='right', va='center', transform=ax.transAxes,
            zorder=4)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def fit_residuals(drpver, redux_path, dapver, analysis_path, daptype, plt, ifu):

    plan_qa_dir = defaults.dap_method_path(daptype, plate=plt, ifudesign=ifu, qa=True,
                                           drpver=drpver, dapver=dapver,
                                           analysis_path=analysis_path)
    if not os.path.isdir(plan_qa_dir):
        os.makedirs(plan_qa_dir)

    rpath = os.path.join(redux_path, str(plt))
    drp_file = os.path.join(rpath, 'stack', 'manga-{0}-{1}-LOGCUBE.fits.gz'.format(plt,ifu))
    image_file = os.path.join(rpath, 'images', '{0}.png'.format(ifu))

    apath = os.path.join(analysis_path, daptype, str(plt), str(ifu))
    maps_file = os.path.join(apath, 'manga-{0}-{1}-MAPS-{2}.fits.gz'.format(plt,ifu,daptype))
    cube_file = os.path.join(apath, 'manga-{0}-{1}-LOGCUBE-{2}.fits.gz'.format(plt,ifu,daptype))

    # Stellar continuum fit qa plot
    ofile_root = os.path.join(plan_qa_dir, 'manga-{0}-{1}-LOGCUBE-{2}-sc-fitqa'.format(
                                plt, ifu, daptype))
    wave, flux, error, model, snrg, rms, frms, rchi2 \
            = get_stellar_continuum_data(maps_file, cube_file)
    ofile = '{0}-lambda.png'.format(ofile_root)
    fom_lambda(plt, ifu, wave, flux, error, model, fit='sc', ofile=ofile)
    ofile = '{0}-growth.png'.format(ofile_root)
    fom_growth(plt, ifu, flux, error, model, snrg, rms, frms, rchi2, fit='sc', ofile=ofile)
    extent, snrg, gmr, rms, frms, rchi2, chi_grw \
            = get_stellar_continuum_fom_maps(drp_file, maps_file)
    ofile = '{0}-maps.png'.format(ofile_root)
    fom_maps(plt, ifu, image_file, snrg, gmr, rms, frms, rchi2, chi_grw, extent=extent, fit='sc',
             ofile=ofile)

    # Emission-line fit qa plot
    ofile_root = os.path.join(plan_qa_dir, 'manga-{0}-{1}-LOGCUBE-{2}-el-fitqa'.format(
                                plt, ifu, daptype))
    wave, flux, error, model, snrg, rms, frms, rchi2 \
            = get_emission_line_data(drp_file, daptype, maps_file, cube_file)
    ofile = '{0}-lambda.png'.format(ofile_root)
    fom_lambda(plt, ifu, wave, flux, error, model, fit='el', ofile=ofile)
    ofile = '{0}-growth.png'.format(ofile_root)
    fom_growth(plt, ifu, flux, error, model, snrg, rms, frms, rchi2, fit='el', ofile=ofile)
    extent, anr, aob, rms, frms, rchi2, chi_grw \
            = get_emission_line_fom_maps(maps_file)
    ofile = '{0}-maps.png'.format(ofile_root)
    fom_maps(plt, ifu, image_file, anr, aob, rms, frms, rchi2, chi_grw, extent=extent, fit='el',
             ofile=ofile)

class FitResiduals(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Construct QA plots showing the fit residuals',
                                    width=width)

        parser.add_argument('plate', type=int, help='plate ID to process')
        parser.add_argument('ifudesign', type=int, help='IFU design to process')

        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument("--redux_path", type=str, help="main DRP output path", default=None)
        parser.add_argument("--analysis_path", type=str, help="main DAP output path", default=None)

        parser.add_argument("--plan_file", type=str, help="parameter file with the MaNGA DAP "
                            "execution plan to use instead of the default" , default=None)

        parser.add_argument('--daptype', type=str, help='DAP processing type', default=None)
        parser.add_argument('--normal_backend', dest='bgagg', action='store_false', default=True)

        return parser

    @staticmethod
    def main(args):
        t = time.perf_counter()

        if args.bgagg:
            pyplot.switch_backend('agg')

        # Set the paths
        redux_path = defaults.drp_redux_path(drpver=args.drpver) \
                        if args.redux_path is None else args.redux_path
        analysis_path = defaults.dap_analysis_path(drpver=args.drpver, dapver=args.dapver) \
                                if args.analysis_path is None else args.analysis_path

        daptypes = []
        if args.daptype is None:
            analysisplan = AnalysisPlanSet.default() if args.plan_file is None \
                            else AnalysisPlanSet.from_par_file(args.plan_file)
            for p in analysisplan:
                bin_method = SpatiallyBinnedSpectra.define_method(p['bin_key'])
                sc_method = StellarContinuumModel.define_method(p['continuum_key'])
                el_method = EmissionLineModel.define_method(p['elfit_key'])
                daptypes += [defaults.dap_method(bin_method['key'],
                                                 sc_method['fitpar']['template_library_key'],
                                                 el_method['continuum_tpl_key'])]
        else:
            daptypes = [args.daptype]

        for daptype in daptypes:
            fit_residuals(args.drpver, redux_path, args.dapver, analysis_path, daptype, args.plate,
                        args.ifudesign)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

