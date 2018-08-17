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
# ax.grid(True, which='major', color='0.8', zorder=1, linestyle='-')
# ax.grid(True, which='minor', color='0.85', zorder=1, linestyle=':')
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


import numpy
from os import remove
import os.path
from time import clock

from mangadap.util.instrument import resample_vector #, log_rebin

from matplotlib import pyplot

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

#    x = numpy.arange(20)
    delt = 100
#    delt = 0
    x = numpy.arange(20)*0.5 + delt
    y = numpy.linspace(0.1, 2.0, 20)
#    y = numpy.square(numpy.sin(50*numpy.linspace(0.1, 1.0, 10)))
    y = y*y - y*y*y + 0.2*y + 3 + 0.2 * numpy.square(numpy.sin(50*y))


    newRange = numpy.array([2,8]) + delt
#    print(x)
#    print(y)

#    X, Y = resample_vector(y, xRange=[x[0],x[-1]], newRange=[4.625,9.375], newpix=20, newLog=False,
#                           flat=False)

#    X, Y = resample_vector(y, xRange=[x[0],x[-1]], newRange=[4,9], newpix=20,
##                           newLog=False,
#                           newLog=True,
#                           flat=False)

#    X, Y = resample_vector(y, xRange=[x[0],x[-1]], newRange=[9, 18.5], newpix=5, newLog=False,
#                           flat=False, ext_value=y[-1])

    Xf, Yf = resample_vector(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                           flat=True)

    X, Y = resample_vector(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                           flat=False)

#    Ym, Xm, scale = log_rebin([x[0],x[-1]], y, log10=True, newRange=newRange, wave_in_ang=True, oversample=6)

#    print(X)
#    print(Y)

    pyplot.scatter(x, y)
    pyplot.plot(x, y)
    pyplot.scatter(X, Y, color='g')
    pyplot.plot(X, Y, color='g')
    pyplot.scatter(Xf, Yf, color='r')
    pyplot.plot(Xf, Yf, color='r')
#    pyplot.scatter(Xm, Ym, color='k')
#    pyplot.plot(Xm, Ym, color='k')
    pyplot.show()

    print('Elapsed time: {0} seconds'.format(clock() - t))



