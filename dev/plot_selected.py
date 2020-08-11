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
# ax = fig.add_axes([0.07, 0.66, 0.3, 0.3], facecolor='0.95')
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
#-----------------------------------------------------------------------
# MATPLOTLIB ERRORBARS
# ax.errorbar(x, y, yerr=ye, xerr=xe, capsize=0, ecolor='0.8', zorder=1)

import os
import time
import numpy

from astropy.io import fits

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from matplotlib import pyplot

from mangadap.drpfits import DRPFitsBitMask

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    select_file = 'manga_d4000_sig_haew_snr_selection.fits'
    spec_file = 'benchmark_spectra_d4000_sig_haew_snr.fits.gz'
    flag_file = 'd4000_sigma_haew_snr_include.db'
    plot_file = 'd4000_sigma_haew_snr_selected.pdf'

    select_hdu = fits.open(select_file)
    spec_hdu = fits.open(spec_file)

    nspec = spec_hdu['FLUX'].shape[0]

    restwave_limits = numpy.array([3575, 7400])

    bm = DRPFitsBitMask()
    flux = numpy.ma.MaskedArray(spec_hdu['FLUX'].data, mask=bm.flagged(spec_hdu['MASK'].data,
                                flag=[ 'DONOTUSE', 'FORESTAR']))

    if not os.path.isfile(flag_file):
        include_flags = numpy.ones(nspec, dtype=int)
        numpy.savetxt(flag_file, numpy.array([ select_hdu['PAR'].data['PLT'].astype(int),
                                               select_hdu['PAR'].data['IFU'].astype(int),
                                               select_hdu['PAR'].data['BIN'].astype(int),
                                               numpy.ones(nspec, dtype=int) ]).T,
                      fmt=[ '%5d', '%5d', '%5d', '%3d' ],
                      header='{0:>3} {1:>5} {2:>5} {3:>3}'.format('PLT', 'IFU', 'BIN', 'SEL'))

    with PdfPages(plot_file) as pdf:

        for i in range(nspec):
#        for i in range(5):
            print('{0}/{1}'.format((i+1),nspec), end='\r')
            lambda_limits = restwave_limits*(1+select_hdu['PAR'].data['NSA_Z'][i])
            indx = numpy.invert(flux.mask[i,:]) & (spec_hdu['WAVE'].data > lambda_limits[0]) \
                        & (spec_hdu['WAVE'].data < lambda_limits[1])
            if numpy.sum(indx) == 0:
                continue
            srt=numpy.argsort(flux[i,indx])
            Df = (flux[i,indx][srt[int(0.9*numpy.sum(indx))]]
                    - flux[i,indx][srt[int(0.1*numpy.sum(indx))]])*3
            flux_limits = [ flux[i,indx][srt[int(0.5*numpy.sum(indx))]] - Df/2, 0]
            flux_limits[1] = flux_limits[0] + Df

            fig = pyplot.figure()

            ax = fig.add_axes([0.15, 0.3, 0.8, 0.6])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10)
            ax.tick_params(which='minor', length=5)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.plot(spec_hdu['WAVE'].data, flux[i,:], zorder=2, color='k', lw=0.5)
            ax.plot(spec_hdu['WAVE'].data, spec_hdu['CONT'].data[i,:], zorder=3, color='C3', lw=0.5)

            ax.set_xlabel(r'$\lambda$')
            ax.set_ylabel(r'$F(\lambda)$')

            ax.text(0.01, 1.03,
                    r'{0}: {1}-{2}-{3}; S/N={4:.1f}; $\sigma$={5:.1f}; D4000={6:.1f}; H$\alpha$ EW={7:.1f}'.format(
                            i, select_hdu['PAR'].data['PLT'][i], select_hdu['PAR'].data['IFU'][i],
                            select_hdu['PAR'].data['BIN'][i], select_hdu['PAR'].data['SNR'][i],
                            select_hdu['PAR'].data['SIGMA'][i], select_hdu['PAR'].data['D4000'][i],
                            select_hdu['PAR'].data['HAEW'][i]),
                    horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

            pdf.savefig(orientation='landscape')
            fig.clear()
            pyplot.close(fig)
        print('{0}/{0}'.format(nspec))

    print('Elapsed time: {0} seconds'.format(time.clock() - t))



