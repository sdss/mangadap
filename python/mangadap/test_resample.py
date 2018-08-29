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
#-----------------------------------------------------------------------
# MATPLOTLIB Logarithmic color bar:
# pcm = ax[0].pcolor(X, Y, Z1,
#                    norm=colors.LogNorm(vmin=Z1.min(), vmax=Z1.max()),
#                    cmap='PuBu_r')
# fig.colorbar(pcm, ax=ax[0], extend='max')
#-----------------------------------------------------------------------
# MATPLOTLIB log-binned 2d map with log-scale axis
# From: https://stackoverflow.com/questions/29175093/creating-a-log-linear-plot-in-matplotlib-using-hist2d
# np.random.seed(1977)
# x, y = np.random.random((2, 1000))
# x = 10**x
# xbins = 10**np.linspace(0, 1, 10)
# ybins = np.linspace(0, 1, 10)
# counts, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
# fig, ax = plt.subplots()
# ax.pcolormesh(xbins, ybins, counts)
# ax.set_xscale('log')
# plt.show()
#-----------------------------------------------------------------------
# PdfPages:
# from matplotlib.backends.backend_pdf import PdfPages
# with PdfPages(ofile) as pdf:
#     for f in ...:
#         ...
#         pdf.savefig(orientation='landscape')
#         fig.clear()
#         pyplot.close()
# return

import os
import time
import numpy

from matplotlib import pyplot

from astropy.io import fits

#from mangadap.util.instrument import resample1d, _pixel_borders
from mangadap.util.instrument import Resample, _pixel_borders
from mangadap.proc.bandpassfilter import passband_integral

import spectres

#-----------------------------------------------------------------------------

def resample_test():
    root = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'])
   
    pltifu = '7815-1901'
    hdu = fits.open(os.path.join(root, pltifu.split('-')[0], 'stack',
                                 'manga-{0}-LOGCUBE.fits.gz'.format(pltifu)))

    drpall_file = os.path.join(root, 'drpall-{0}.fits'.format(os.environ['MANGADRP_VER']))
    drpall = fits.open(drpall_file)[1].data

    indx = drpall['PLATEIFU'] == pltifu
    z = drpall['NSA_Z'][indx][0]
    print(z)

    old_wave = hdu['WAVE'].data
    old_flux = hdu['FLUX'].data[:,10,10]
    indx = (old_wave > old_wave[0]/(1+z)) & (old_wave < old_wave[-2]/(1+z))

    pyplot.plot(old_wave/(1+z), old_flux)

    t = time.clock()
    new_flux_spectres = spectres.spectres(old_wave[indx], old_wave/(1+z), old_flux)
    print('SpectRes Time: ', time.clock()-t)

    pyplot.plot(old_wave[indx], new_flux_spectres)

    t = time.clock()
    borders = _pixel_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size, log=True)[0]
    _p = numpy.repeat(borders, 2)[1:-1].reshape(-1,2)
    new_flux_brute = passband_integral(old_wave/(1+z), old_flux, passband=_p,
                                       log=True)/(_p[:,1]-_p[:,0])
    print('Brute Force Time: ', time.clock()-t)

    pyplot.plot(old_wave, new_flux_brute)

    t = time.clock()
    r = Resample(old_flux, x=old_wave/(1+z), newRange=[old_wave[0], old_wave[-1]], inLog=True,
                 newLog=True)
    print('Resample Time: ', time.clock()-t)

    pyplot.plot(r.outx, r.outy)

    print('Mean diff:')
    print('    spectres - brute    = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_spectres-new_flux_brute[indx]))))
    print('    spectres - resample = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_spectres-r.outy[indx]))))
    print('    brute - resample    = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_brute-r.outy))))

    pyplot.show()



#-----------------------------------------------------------------------------

if __name__ == '__main__':
    resample_test()



