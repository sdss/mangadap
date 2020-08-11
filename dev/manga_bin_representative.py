#!/usr/bin/env python3

# Useful stuff:
#-----------------------------------------------------------------------
# ASTROPY CONSTANTS
# from astropy import constants
# constants.c.to('km/s').value
#-----------------------------------------------------------------------
# MATPLOTLIB TICKMARKS: ORDERING WILL MATTER TO OUTPUT
# from matplotlib.ticker import NullFormatter
# ax.xaxis.set_major_formatter(NullFormatter())
# ax.minorticks_on()
# ax.tick_params(which='major', length=6)
# ax.tick_params(which='minor', length=3)
#
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# majorLocator   = MultipleLocator(20)
# majorFormatter = FormatStrFormatter('%d')
# minorLocator   = MultipleLocator(5)
# ax.xaxis.set_major_locator(majorLocator)
# ax.xaxis.set_major_formatter(majorFormatter)
# ax.xaxis.set_minor_locator(minorLocator)
#-----------------------------------------------------------------------
# MATPLOTLIB PRINT TO A FILE
# fig = pyplot.figure(1)
# fig.canvas.print_figure('out.pdf')
# fig.canvas.print_figure('test.pdf', bbox_inches='tight')
#-----------------------------------------------------------------------
# MATPLOTLIB LOG SCALE
# pyplot.xscale('log')
#-----------------------------------------------------------------------
# MATPLOTLIB TEXT IN PLOT VIEWPORT COORDINATES
# ax.text(-0.25, 0.50, r'$\Delta V$ (MIUSCAT-MILES)', horizontalalignment='center',
#         verticalalignment='center', transform=ax.transAxes, rotation='vertical')
#-----------------------------------------------------------------------
# MATPLOTLIB SET FIGURE ASPECT RATIO
# w,h = pyplot.figaspect(1)
# fig = pyplot.figure(figsize=(1.5*w,1.5*h))
#-----------------------------------------------------------------------
# MATPLOTLIB CREATE AXIS WITH VIEWPORT COORDINATES
# ax = pyplot.axes([ left, bottom+i*height, width, height ])
#   OR IN A FIGURE (with a gray background):
# ax = fig.add_axes([0.07, 0.66, 0.3, 0.3], axisbg='0.95')
#-----------------------------------------------------------------------
# MATPLOTLIB ADD GRID TO PANEL
# ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
# ax.grid(True, which='minor', color='0.85', zorder=0, linestyle=':')
#-----------------------------------------------------------------------
# MATPLOTLIB CREATE PLOT PANEL GRID
# from matplotlib import gridspec
# gs = gridspec.GridSpec(100,100)
# ax = pyplot.subplot(gs[0:30,4:34])
#-----------------------------------------------------------------------
# MATPLOTLIB CHANGE FONT (OR OTHER PROPERTIES)
# from matplotlib import rc
# font = { 'size' : 8 }
# rc('font', **font)
#-----------------------------------------------------------------------
# MATPLOTLIB CREATE SCATTER PLOT WITH COLOR BAR
# import matplotlib.pyplot as plt
# cm = plt.cm.get_cmap('RdYlBu')
# xy = range(20)
# z = xy
# sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=20, s=35, cmap=cm)
# plt.colorbar(sc)
# plt.show()
#-----------------------------------------------------------------------
# MATPLOTLIB SHADED BOX
# box_sig_x = [ sig_lim[0], sigi_x_1901, sigi_x_1901, sig_lim[0] ]
# box_sig_y = [ sig_lim[0], sig_lim[0], sigi_y_1901, sigi_y_1901 ]
# ax.fill_between( box_sig_x, box_sig_y, color='blue', alpha=0.1, zorder=1)
#-----------------------------------------------------------------------
# MATPLOTLIB OPPOSITE AXIS LABELS, LIMITS
# axt = ax.twinx()
# axt.set_xlim(alim)
# axt.set_ylim([1,45])
# axt.minorticks_on()
# axt.tick_params(which='major', length=10)
# axt.tick_params(which='minor', length=5)
# axt.xaxis.set_major_locator(MultipleLocator(0.05))
# axt.xaxis.set_minor_locator(MultipleLocator(0.01))
# axt.yaxis.set_major_locator(MultipleLocator(20))
# axt.yaxis.set_major_formatter(FormatStrFormatter('%d'))
# axt.yaxis.set_minor_locator(MultipleLocator(10))
# axt.xaxis.set_major_formatter(NullFormatter())
# axt.text(1.2, 0.5, r'$P\ (a)$', horizontalalignment='center', verticalalignment='center',
#             transform=axt.transAxes, rotation='vertical')
#-----------------------------------------------------------------------
# PARSING A STRING MATHEMATICAL EXPRESSION
# import parser
# import numpy
# t = 2
# st = parser.expr('t**3+4')
# code = st.compile()
# eval(code)
# t = numpy.arange(3)
# eval(code)
#-----------------------------------------------------------------------
# MATPLOTLIB DISCRETE COLOR BAR
# from matplotlib import cm
# cmap = cm.get_cmap('inferno', 6)

import os
import time
import numpy
import glob

from astropy.io import fits
import astropy.constants

#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot, rc, ticker, colors

from mangadap.dapfits import DAPQualityBitMask, DAPMapsBitMask
from mangadap.config.defaults import default_dap_method_path

#-----------------------------------------------------------------------------

def twod_bin(x, y, x0, dx, nx, y0, dy, ny, z=None):
    _z = numpy.ones(x.shape, dtype=float) if z is None else z

    bini = ((x-x0)/dx).astype(int)
    bini[ (bini < 0) | (bini >= nx) ] = -1
    binj = ((y-y0)/dy).astype(int)
    binj[ (binj < 0) | (binj >= ny) ] = -1

    n = numpy.zeros((nx,ny),dtype=int)
#    m = numpy.zeros((nx,ny),dtype=float)
#    d = numpy.zeros((nx,ny),dtype=float)
    for i in range(nx):
        iindx = bini == i
        for j in range(ny):
            jindx = iindx & (binj == j)
            if numpy.sum(jindx) == 0:
                continue
            n[i,j] = numpy.sum(jindx)
#            m[i,j] = numpy.median(_z[jindx])
#            d[i,j] = median_absolute_deviation(_z[jindx])

    return n if z is None else (n,m,d)


def d4000_sigma_haew_dist(fits_file, n, d4000_lim, sigma_lim, haew_lim, min_sn=None):

    # Channels with the H-alpha data and the D4000 data
    hac = 18
    d4c = 43

    # Set the method and get the root directory
    method = 'SPX-GAU-MILESHC'
    odir = default_dap_method_path(method)
    print('Output directory: {0}'.format(odir))

    # Find the intermediate files with the stellar-continuum models
    srchstr = os.path.join(odir, '*', '*', 'manga-*-MAPS-{0}.fits.gz'.format(method))
    print('Search string: {0}'.format(srchstr))
    maps_files = glob.glob(srchstr)
    nfiles = len(maps_files)
    print('Found {0} MAPS files.'.format(nfiles))

    d4000_bins = numpy.linspace(*d4000_lim,n+1)
    sigma_bins = numpy.logspace(*numpy.log10(sigma_lim),n+1)
    haew_bins = numpy.logspace(*numpy.log10(haew_lim),n+1)

    dl = d4000_bins[0]
    dd = d4000_bins[1]- d4000_bins[0]
    sl = numpy.log10(sigma_bins[0])
    ds = numpy.mean(numpy.diff(numpy.log10(sigma_bins)))
    hl = numpy.log10(haew_bins[0])
    dh = numpy.mean(numpy.diff(numpy.log10(haew_bins)))

    d4000_vs, vs_sigma = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                        sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2)
    d4000_vs, vs_haew = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                       haew_bins[:-1]+dh*haew_bins[:-1]*numpy.log(10)/2)
    sigma_vs, vs_haew = numpy.meshgrid(sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2,
                                        haew_bins[:-1]+ds*haew_bins[:-1]*numpy.log(10)/2)

    bin_d4000_vs_sigma = numpy.zeros((n,n), dtype=int)
    bin_d4000_vs_haew = numpy.zeros((n,n), dtype=int)
    bin_sigma_vs_haew = numpy.zeros((n,n), dtype=int)

    qual_bm = DAPQualityBitMask()
    maps_bm = DAPMapsBitMask()

    for fi in range(nfiles):
#    for fi in range(200):
        print('{0}/{1}'.format(fi+1,nfiles), end='\r')
        maps_file = maps_files[fi]
        print(maps_file)

        # Get the bin IDs, mapped to the correct spaxels
        hdu = fits.open(maps_file)

        # Check if the file was CRITICAL
        if qual_bm.flagged(hdu['PRIMARY'].header['DAPQUAL'], flag='DRPCRIT'):
            print('DRPCRIT')
            print(maps_file)
            hdu.close()
            del hdu
            continue

        plt, ifu = tuple([ int(x) for x in hdu[0].header['PLATEIFU'].split('-')])
        nsa_z = hdu['PRIMARY'].header['SCINPVEL'] / astropy.constants.c.to('km/s').value
        binid_map = hdu['BINID'].data.copy()

        spatial_shape = binid_map.shape
        binid, bin_indx = numpy.unique(binid_map.ravel(), return_index=True)
        binid = binid[1:]
        bin_indx = bin_indx[1:]

        # Get the D4000, H-alpha EW, and sigma values
#        snr = hdu['SPX_SNR'].data.ravel()[bin_indx]
#        snr_indx = snr > 1
#        good_snr = snr[snr_indx]
#        binid_snr = binid[snr_indx]

        mask = (hdu['EMLINE_GEW'].data[hac,:,:].ravel()[bin_indx] > 0) \
               & (hdu['SPECINDEX'].data[d4c,:,:].ravel()[bin_indx] > 0) \
               & (hdu['STELLAR_SIGMA'].data.ravel()[bin_indx] > 0)
        if min_sn is not None:
            mask &= (hdu['SPX_SNR'].data.ravel()[bin_indx] > min_sn)

        haew = hdu['EMLINE_GEW'].data[hac,:,:].ravel()[bin_indx][mask] #[snr_indx]
        d4000 = hdu['SPECINDEX'].data[d4c,:,:].ravel()[bin_indx][mask] #[snr_indx]
        sigma = hdu['STELLAR_SIGMA'].data.ravel()[bin_indx][mask] #[snr_indx]

        hdu.close()
        del hdu

#        bin_d4000_vs_sigma += twod_bin(numpy.log10(d4000), numpy.log10(sigma),
        bin_d4000_vs_sigma += twod_bin(d4000, numpy.log10(sigma),
                                       dl, dd, n, sl, ds, n)
#        bin_d4000_vs_haew += twod_bin(numpy.log10(d4000), numpy.log10(haew),
        bin_d4000_vs_haew += twod_bin(d4000, numpy.log10(haew),
                                      dl, dd, n, hl, dh, n)
        bin_sigma_vs_haew += twod_bin(numpy.log10(sigma), numpy.log10(haew),
                                      sl, ds, n, hl, dh, n)

    fits.HDUList([ fits.PrimaryHDU(),
                   fits.ImageHDU(data=bin_d4000_vs_sigma, name='D4000SIGMA'),
                   fits.ImageHDU(data=bin_d4000_vs_haew, name='D4000HAEW'),
                   fits.ImageHDU(data=bin_sigma_vs_haew, name='SIGMAHAEW')
                  ]).writeto(fits_file, overwrite=True)


def plot_dist_3d(fits_file, n, d4000_lim, sigma_lim, haew_lim):

    # Set the bin limits (not read from file!!)
    d4000_bins = numpy.linspace(*d4000_lim,n+1)
    sigma_bins = numpy.logspace(*numpy.log10(sigma_lim),n+1)
    haew_bins = numpy.logspace(*numpy.log10(haew_lim),n+1)

    dl = d4000_bins[0]
    dd = d4000_bins[1]- d4000_bins[0]
    sl = numpy.log10(sigma_bins[0])
    ds = numpy.mean(numpy.diff(numpy.log10(sigma_bins)))
    hl = numpy.log10(haew_bins[0])
    dh = numpy.mean(numpy.diff(numpy.log10(haew_bins)))

    d4000_vs, vs_sigma = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                        sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2)
    d4000_vs, vs_haew = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                       haew_bins[:-1]+dh*haew_bins[:-1]*numpy.log(10)/2)
    sigma_vs, vs_haew = numpy.meshgrid(sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2,
                                        haew_bins[:-1]+ds*haew_bins[:-1]*numpy.log(10)/2)

    hdu = fits.open(fits_file)

    bin_d4000_vs_sigma = hdu['D4000SIGMA'].data.copy()
    bin_d4000_vs_haew = hdu['D4000HAEW'].data.copy()
    bin_sigma_vs_haew = hdu['SIGMAHAEW'].data.copy()

    font = { 'size' : 12 }
    rc('font', **font)

    n_lim = numpy.amax([ numpy.amax(bin_d4000_vs_sigma), numpy.amax(bin_d4000_vs_haew),
                  numpy.amax(bin_sigma_vs_haew) ])
    print(n_lim)
    levels = [n_lim/32, n_lim/16, n_lim/8, n_lim/4, n_lim/2]

    _bin_d4000_vs_sigma = numpy.ma.MaskedArray(bin_d4000_vs_sigma.T,
                                               mask=numpy.invert(bin_d4000_vs_sigma.T>0))
    _bin_d4000_vs_haew = numpy.ma.MaskedArray(bin_d4000_vs_haew.T,
                                              mask=numpy.invert(bin_d4000_vs_haew.T>0))
    _bin_sigma_vs_haew = numpy.ma.MaskedArray(bin_sigma_vs_haew.T,
                                               mask=numpy.invert(bin_sigma_vs_haew.T>0))
    
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    cax = fig.add_axes([0.66, 0.66, 0.2, 0.02])

    ax = fig.add_axes([0.15, 0.56, 0.4, 0.4], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(d4000_lim)
    ax.set_ylim(sigma_lim)
    ax.set_yscale('log')
    ax.grid(True, which='major', color='0.85', zorder=1, linestyle='-')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
#    ax.set_xscale('log')
#    img = ax.imshow(bin_d4000_vs_sigma.T, origin='lower', interpolation='nearest',
#              extent=extent_d4000_vs_sigma, aspect='auto', zorder=2, vmin=0, vmax=n_lim)
    img = ax.pcolormesh(d4000_bins, sigma_bins, _bin_d4000_vs_sigma,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(d4000_vs, vs_sigma, _bin_d4000_vs_sigma, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    pyplot.colorbar(img, cax=cax, orientation='horizontal') #, ticks=ticker.MultipleLocator(200))
    ax.text(-0.2, 0.5, r'$\sigma_{\rm obs}$ (km/s)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')
    cax.text(0.5, -3, r'$N_{\rm spaxel}$', horizontalalignment='center',
             verticalalignment='center', transform=cax.transAxes)

    cax.text(0.5, 8, r'MPL-6', horizontalalignment='center',
             verticalalignment='center', transform=cax.transAxes, fontsize=20)
    cax.text(0.5, 6, r'S/N > 10', horizontalalignment='center',
             verticalalignment='center', transform=cax.transAxes, fontsize=20)

    ax = fig.add_axes([0.15, 0.15, 0.4, 0.4], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(d4000_lim)
    ax.set_ylim(haew_lim)
    ax.set_yscale('log')
#    ax.set_xscale('log')
#    ax.imshow(bin_d4000_vs_haew.T, origin='lower', interpolation='nearest',
#              extent=extent_d4000_vs_haew, aspect='auto', zorder=2, vmin=0, vmax=n_lim)
    img = ax.pcolormesh(d4000_bins, haew_bins, _bin_d4000_vs_haew,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(d4000_vs, vs_haew, _bin_d4000_vs_haew, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    ax.text(-0.2, 0.5, r'H$\alpha$ EW ($\AA$)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')
    ax.text(0.5, -0.15, r'D4000', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)

    ax = fig.add_axes([0.56, 0.15, 0.4, 0.4], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(sigma_lim)
    ax.set_ylim(haew_lim)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    ax.imshow(bin_d4000_vs_haew.T, origin='lower', interpolation='nearest',
#              extent=extent_sigma_vs_haew, aspect='auto', zorder=2, vmin=0, vmax=n_lim)
    img = ax.pcolormesh(sigma_bins, haew_bins, _bin_sigma_vs_haew,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(sigma_vs, vs_haew, _bin_sigma_vs_haew, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    ax.text(0.5, -0.15, r'$\sigma_{\rm obs}$ (km/s)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)

    plot_file = 'manga_representative.pdf'
    plot_file = None
    if plot_file is None:
        pyplot.show()
    else:
        fig.canvas.print_figure(plot_file, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)
    

    print('DONE                     ')


# Include binning with S/N
def d4000_sigma_haew_snr_dist(fits_file, n, d4000_lim, sigma_lim, haew_lim, snr_lim):

    # Channels with the H-alpha data and the D4000 data
    hac = 18
    d4c = 43

    # Set the method and get the root directory
    method = 'SPX-GAU-MILESHC'
    odir = default_dap_method_path(method)
    print('Output directory: {0}'.format(odir))

    # Find the intermediate files with the stellar-continuum models
    srchstr = os.path.join(odir, '*', '*', 'manga-*-MAPS-{0}.fits.gz'.format(method))
    print('Search string: {0}'.format(srchstr))
    maps_files = glob.glob(srchstr)
    nfiles = len(maps_files)
    print('Found {0} MAPS files.'.format(nfiles))

    d4000_bins = numpy.linspace(*d4000_lim,n+1)
    sigma_bins = numpy.logspace(*numpy.log10(sigma_lim),n+1)
    haew_bins = numpy.logspace(*numpy.log10(haew_lim),n+1)
    snr_bins = numpy.logspace(*numpy.log10(snr_lim),n+1)

    dl = d4000_bins[0]
    dd = d4000_bins[1]- d4000_bins[0]
    sl = numpy.log10(sigma_bins[0])
    ds = numpy.mean(numpy.diff(numpy.log10(sigma_bins)))
    hl = numpy.log10(haew_bins[0])
    dh = numpy.mean(numpy.diff(numpy.log10(haew_bins)))
    nl = numpy.log10(snr_bins[0])
    dn = numpy.mean(numpy.diff(numpy.log10(snr_bins)))

    d4000_vs, vs_sigma = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                        sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2)
    d4000_vs, vs_haew = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                       haew_bins[:-1]+dh*haew_bins[:-1]*numpy.log(10)/2)
    d4000_vs, vs_snr = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                      snr_bins[:-1]+dn*snr_bins[:-1]*numpy.log(10)/2)
    sigma_vs, vs_haew = numpy.meshgrid(sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2,
                                        haew_bins[:-1]+ds*haew_bins[:-1]*numpy.log(10)/2)
    sigma_vs, vs_snr = numpy.meshgrid(sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2,
                                      snr_bins[:-1]+dn*snr_bins[:-1]*numpy.log(10)/2)
    haew_vs, vs_snr = numpy.meshgrid(haew_bins[:-1]+dh*haew_bins[:-1]*numpy.log(10)/2,
                                     snr_bins[:-1]+dn*snr_bins[:-1]*numpy.log(10)/2)

    bin_d4000_vs_sigma = numpy.zeros((n,n), dtype=int)
    bin_d4000_vs_haew = numpy.zeros((n,n), dtype=int)
    bin_d4000_vs_snr = numpy.zeros((n,n), dtype=int)
    bin_sigma_vs_haew = numpy.zeros((n,n), dtype=int)
    bin_sigma_vs_snr = numpy.zeros((n,n), dtype=int)
    bin_haew_vs_snr = numpy.zeros((n,n), dtype=int)

    qual_bm = DAPQualityBitMask()
    maps_bm = DAPMapsBitMask()

    for fi in range(nfiles):
#    for fi in range(20):
        print('{0}/{1}'.format(fi+1,nfiles), end='\r')
        maps_file = maps_files[fi]
        print(maps_file)

        # Get the bin IDs, mapped to the correct spaxels
        hdu = fits.open(maps_file)

        # Check if the file was CRITICAL
        if qual_bm.flagged(hdu['PRIMARY'].header['DAPQUAL'], flag='DRPCRIT'):
            print('DRPCRIT')
            print(maps_file)
            hdu.close()
            del hdu
            continue

        plt, ifu = tuple([ int(x) for x in hdu[0].header['PLATEIFU'].split('-')])
        nsa_z = hdu['PRIMARY'].header['SCINPVEL'] / astropy.constants.c.to('km/s').value
        binid_map = hdu['BINID'].data[0,:,:].copy()

        spatial_shape = binid_map.shape
        binid, bin_indx = numpy.unique(binid_map.ravel(), return_index=True)
        binid = binid[1:]
        bin_indx = bin_indx[1:]

        # Get the D4000, H-alpha EW, and sigma values
#        snr = hdu['SPX_SNR'].data.ravel()[bin_indx]
#        snr_indx = snr > 1
#        good_snr = snr[snr_indx]
#        binid_snr = binid[snr_indx]

        mask = (hdu['EMLINE_GEW'].data[hac,:,:].ravel()[bin_indx] > 0) \
               & (hdu['SPECINDEX'].data[d4c,:,:].ravel()[bin_indx] > 0) \
               & (hdu['STELLAR_SIGMA'].data.ravel()[bin_indx] > 0) \
               & (hdu['SPX_SNR'].data.ravel()[bin_indx] > 0)

        haew = hdu['EMLINE_GEW'].data[hac,:,:].ravel()[bin_indx][mask]
        d4000 = hdu['SPECINDEX'].data[d4c,:,:].ravel()[bin_indx][mask]
        sigma = hdu['STELLAR_SIGMA'].data.ravel()[bin_indx][mask]
        snr = hdu['SPX_SNR'].data.ravel()[bin_indx][mask]

        hdu.close()
        del hdu

        bin_d4000_vs_sigma += twod_bin(d4000, numpy.log10(sigma), dl, dd, n, sl, ds, n)
        bin_d4000_vs_haew += twod_bin(d4000, numpy.log10(haew), dl, dd, n, hl, dh, n)
        bin_d4000_vs_snr += twod_bin(d4000, numpy.log10(snr), dl, dd, n, nl, dn, n)
        bin_sigma_vs_haew += twod_bin(numpy.log10(sigma), numpy.log10(haew), sl, ds, n, hl, dh, n)
        bin_sigma_vs_snr += twod_bin(numpy.log10(sigma), numpy.log10(snr), sl, ds, n, nl, dn, n)
        bin_haew_vs_snr += twod_bin(numpy.log10(haew), numpy.log10(snr), hl, dh, n, nl, dn, n)

    fits.HDUList([ fits.PrimaryHDU(),
                   fits.ImageHDU(data=bin_d4000_vs_sigma, name='D4000SIGMA'),
                   fits.ImageHDU(data=bin_d4000_vs_haew, name='D4000HAEW'),
                   fits.ImageHDU(data=bin_d4000_vs_snr, name='D4000SNR'),
                   fits.ImageHDU(data=bin_sigma_vs_haew, name='SIGMAHAEW'),
                   fits.ImageHDU(data=bin_sigma_vs_snr, name='SIGMASNR'),
                   fits.ImageHDU(data=bin_haew_vs_snr, name='HAEWSNR')
                  ]).writeto(fits_file, overwrite=True)


def plot_dist_4d(fits_file, n, d4000_lim, sigma_lim, haew_lim, snr_lim):

    # Set the bin limits (not read from file!!)
    d4000_bins = numpy.linspace(*d4000_lim,n+1)
    sigma_bins = numpy.logspace(*numpy.log10(sigma_lim),n+1)
    haew_bins = numpy.logspace(*numpy.log10(haew_lim),n+1)
    snr_bins = numpy.logspace(*numpy.log10(snr_lim),n+1)

    dl = d4000_bins[0]
    dd = d4000_bins[1]- d4000_bins[0]
    sl = numpy.log10(sigma_bins[0])
    ds = numpy.mean(numpy.diff(numpy.log10(sigma_bins)))
    hl = numpy.log10(haew_bins[0])
    dh = numpy.mean(numpy.diff(numpy.log10(haew_bins)))
    nl = numpy.log10(snr_bins[0])
    dn = numpy.mean(numpy.diff(numpy.log10(snr_bins)))

    d4000_vs, vs_sigma = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                        sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2)
    d4000_vs, vs_haew = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                       haew_bins[:-1]+dh*haew_bins[:-1]*numpy.log(10)/2)
    d4000_vs, vs_snr = numpy.meshgrid(d4000_bins[:-1]+dd/2,
                                      snr_bins[:-1]+dn*snr_bins[:-1]*numpy.log(10)/2)
    sigma_vs, vs_haew = numpy.meshgrid(sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2,
                                        haew_bins[:-1]+ds*haew_bins[:-1]*numpy.log(10)/2)
    sigma_vs, vs_snr = numpy.meshgrid(sigma_bins[:-1]+ds*sigma_bins[:-1]*numpy.log(10)/2,
                                      snr_bins[:-1]+dn*snr_bins[:-1]*numpy.log(10)/2)
    haew_vs, vs_snr = numpy.meshgrid(haew_bins[:-1]+dh*haew_bins[:-1]*numpy.log(10)/2,
                                     snr_bins[:-1]+dn*snr_bins[:-1]*numpy.log(10)/2)

    hdu = fits.open(fits_file)

    bin_d4000_vs_sigma = hdu['D4000SIGMA'].data.copy()
    bin_d4000_vs_haew = hdu['D4000HAEW'].data.copy()
    bin_d4000_vs_snr = hdu['D4000SNR'].data.copy()
    bin_sigma_vs_haew = hdu['SIGMAHAEW'].data.copy()
    bin_sigma_vs_snr = hdu['SIGMASNR'].data.copy()
    bin_haew_vs_snr = hdu['HAEWSNR'].data.copy()

    n_lim = numpy.amax([ numpy.amax(bin_d4000_vs_sigma),
                         numpy.amax(bin_d4000_vs_haew),
                         numpy.amax(bin_d4000_vs_snr),
                         numpy.amax(bin_sigma_vs_haew),
                         numpy.amax(bin_sigma_vs_snr),
                         numpy.amax(bin_haew_vs_snr) ])
    print(n_lim)
    levels = [n_lim/32, n_lim/16, n_lim/8, n_lim/4, n_lim/2]

    _bin_d4000_vs_sigma = numpy.ma.MaskedArray(bin_d4000_vs_sigma.T,
                                               mask=numpy.invert(bin_d4000_vs_sigma.T>0))
    _bin_d4000_vs_haew = numpy.ma.MaskedArray(bin_d4000_vs_haew.T,
                                              mask=numpy.invert(bin_d4000_vs_haew.T>0))
    _bin_d4000_vs_snr = numpy.ma.MaskedArray(bin_d4000_vs_snr.T,
                                             mask=numpy.invert(bin_d4000_vs_snr.T>0))
    _bin_sigma_vs_haew = numpy.ma.MaskedArray(bin_sigma_vs_haew.T,
                                               mask=numpy.invert(bin_sigma_vs_haew.T>0))
    _bin_sigma_vs_snr = numpy.ma.MaskedArray(bin_sigma_vs_snr.T,
                                             mask=numpy.invert(bin_sigma_vs_snr.T>0))
    _bin_haew_vs_snr = numpy.ma.MaskedArray(bin_haew_vs_snr.T,
                                             mask=numpy.invert(bin_haew_vs_snr.T>0))
    
    font = { 'size' : 12 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    cax = fig.add_axes([0.59, 0.79, 0.2, 0.02])

    ax = fig.add_axes([0.11, 0.68, 0.28, 0.28], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(d4000_lim)
    ax.set_ylim(sigma_lim)
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    img = ax.pcolormesh(d4000_bins, sigma_bins, _bin_d4000_vs_sigma,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(d4000_vs, vs_sigma, _bin_d4000_vs_sigma, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    pyplot.colorbar(img, cax=cax, orientation='horizontal') #, ticks=ticker.MultipleLocator(200))
    ax.text(-0.25, 0.5, r'$\sigma_{\rm obs}$ (km/s)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')
    cax.text(0.5, -3, r'$N_{\rm spaxel}$', horizontalalignment='center',
             verticalalignment='center', transform=cax.transAxes)

    cax.text(0.5, 4, r'MPL-6', horizontalalignment='center',
             verticalalignment='center', transform=cax.transAxes, fontsize=20)

    ax = fig.add_axes([0.11, 0.39, 0.28, 0.28], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(d4000_lim)
    ax.set_ylim(haew_lim)
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    img = ax.pcolormesh(d4000_bins, haew_bins, _bin_d4000_vs_haew,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(d4000_vs, vs_haew, _bin_d4000_vs_haew, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    ax.text(-0.25, 0.5, r'H$\alpha$ EW ($\AA$)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')

    ax = fig.add_axes([0.11, 0.10, 0.28, 0.28], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(d4000_lim)
    ax.set_ylim(snr_lim)
    ax.set_yscale('log')
    img = ax.pcolormesh(d4000_bins, snr_bins, _bin_d4000_vs_snr,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(d4000_vs, vs_snr, _bin_d4000_vs_snr, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    ax.text(-0.25, 0.5, r'S/N$_g$', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')
    ax.text(0.5, -0.2, r'D4000', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)

    ax = fig.add_axes([0.40, 0.39, 0.28, 0.28], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(sigma_lim)
    ax.set_ylim(haew_lim)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    img = ax.pcolormesh(sigma_bins, haew_bins, _bin_sigma_vs_haew,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(sigma_vs, vs_haew, _bin_sigma_vs_haew, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))

    ax = fig.add_axes([0.40, 0.10, 0.28, 0.28], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(sigma_lim)
    ax.set_ylim(snr_lim)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    img = ax.pcolormesh(sigma_bins, snr_bins, _bin_sigma_vs_snr,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(sigma_vs, vs_snr, _bin_sigma_vs_snr, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    ax.text(0.5, -0.20, r'$\sigma_{\rm obs}$ (km/s)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)

    ax = fig.add_axes([0.69, 0.10, 0.28, 0.28], facecolor='0.95')
    ax.minorticks_on()
    ax.set_xlim(haew_lim)
    ax.set_ylim(snr_lim)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    img = ax.pcolormesh(haew_bins, snr_bins, _bin_haew_vs_snr,
                        zorder=3, norm=colors.LogNorm(vmin=1, vmax=n_lim), lw=0, rasterized=True)
    ax.contour(haew_vs, vs_snr, _bin_haew_vs_snr, zorder=4,
               levels=levels, linewidths=0.5, norm=colors.LogNorm(vmin=1/3, vmax=n_lim/3))
    ax.text(0.5, -0.20, r'H$\alpha$ EW ($\AA$)', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)

    plot_file = 'manga_representative_4d.pdf'
#    plot_file = None
    if plot_file is None:
        pyplot.show()
    else:
        fig.canvas.print_figure(plot_file, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)
    

    print('DONE                     ')


#    fits.HDUList([ fits.PrimaryHDU(),
#                   fits.BinTableHDU.from_columns([ fits.Column(name=n,
#                                                               format=rec_to_fits_type(tbl[n]),
#                                                               array=tbl[n])
#                                                       for n in tbl.dtype.names ],
#                                                  name='PAR')
#                  ]).writeto('crap.fits.gz', overwrite=True)
#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    # No S/N limit
    n=100
    d4000_lim=[0.7,3.5]
    sigma_lim=[20,600]
    haew_lim=[0.03,200]
    snr_lim=[1,200]

    fits_file='SPX-bin_snr.fits.gz'
#    overwrite=True
    overwrite=False
    if not os.path.isfile(fits_file) or overwrite:
        d4000_sigma_haew_snr_dist(fits_file, n, d4000_lim, sigma_lim, haew_lim, snr_lim)
    plot_dist_4d(fits_file, n, d4000_lim, sigma_lim, haew_lim, snr_lim)

    exit()

    # No S/N limit
    n=100
    d4000_lim=[0.7,3.5]
    sigma_lim=[20,600]
    haew_lim=[0.03,200]

    fits_file='SPX-snr_gt_00.fits.gz'
#    overwrite=True
    overwrite=False
    if not os.path.isfile(fits_file) or overwrite:
        d4000_sigma_haew_dist(fits_file, n, d4000_lim, sigma_lim, haew_lim)
    plot_dist_3d(fits_file, n, d4000_lim, sigma_lim, haew_lim)

    # Set S/N > 10 limit, but keep the same volume dimensions
    n=100
    d4000_lim=[0.7,3.5]
    sigma_lim=[20,600]
    haew_lim=[0.03,200]

    fits_file='SPX-snr_gt_10.fits.gz'
#    overwrite=True
    overwrite=False
    if not os.path.isfile(fits_file) or overwrite:
        d4000_sigma_haew_dist(fits_file, n, d4000_lim, sigma_lim, haew_lim, min_sn=10)
    plot_dist_3d(fits_file, n, d4000_lim, sigma_lim, haew_lim)

    # Repeat S/N > 10 limit, but with smaller d4000 dimension
    n=100
    d4000_lim=[0.8,2.6]
    sigma_lim=[20,600]
    haew_lim=[0.03,200]

    fits_file='SPX-snr_gt_10_v2.fits.gz'
#    overwrite=True
    overwrite=False
    if not os.path.isfile(fits_file) or overwrite:
        d4000_sigma_haew_dist(fits_file, n, d4000_lim, sigma_lim, haew_lim, mins_sn=10)
    plot_dist_3d(fits_file, n, d4000_lim, sigma_lim, haew_lim)
    print('Elapsed time: {0} seconds'.format(time.clock() - t))


