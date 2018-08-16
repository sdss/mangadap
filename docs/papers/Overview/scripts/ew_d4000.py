#!/usr/bin/env python3

import os
import time
import numpy

from matplotlib import pyplot, ticker, rc, colors

from astropy.io import fits
import astropy.constants

from mangadap.proc.util import growth_lim
from mangadap.util.bitmask import BitMask

#-----------------------------------------------------------------------------
def init_ax(fig, pos):
    ax = fig.add_axes(pos, facecolor='0.95')
    ax.grid(True, which='major', color='0.9', zorder=0, linestyle='-')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    return ax


def ew_d4000(drpall, dapall):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == 'HYB10-GAU-MILESHC')
    indx &= (drpall['nsa_sersic_mass'][dapall['DRPALLINDX']] > 0)

#    ha_indx = indx & (dapall['EMLINE_GSB_1RE'][:,18] > -999) \
#                    & (dapall['SPECINDEX_1RE'][:,44] > -999)
    ha_indx = indx & (dapall['EMLINE_GEW_1RE'][:,18] > -999) \
                    & (dapall['SPECINDEX_1RE'][:,44] > -999)
    hd_indx = indx & (dapall['SPECINDEX_1RE'][:,21] > -999) \
                    & (dapall['SPECINDEX_1RE'][:,44] > -999)

    mass_lim = [1e9, 7e11]
    d4_lim = [0.9, 2.3]
#    ha_lim = [0.005, 800]
    ha_lim = [2e-2, 400]
    hd_lim = [-4, 9]

    font = { 'size' : 12 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.3, 0.56, 0.4, 0.4])
    ax.set_xlim(d4_lim)
    ax.set_ylim(ha_lim)
    ax.tick_params(axis='x', which='both', direction='in')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())

    sc = ax.scatter(dapall['SPECINDEX_1RE'][ha_indx,44], dapall['EMLINE_GEW_1RE'][ha_indx,18],
                    lw=0, marker='.', s=10, zorder=2,
                    c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][ha_indx]],
                    norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
                    cmap='viridis')

    ax.text(-0.25, 0.5, r'$\langle {\rm EW}_{{\rm H}\alpha}\rangle_e$ [$\AA$, in emission]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax = init_ax(fig, [0.3, 0.15, 0.4, 0.4])
    ax.set_xlim(d4_lim)
    ax.set_ylim(hd_lim)

    sc = ax.scatter(dapall['SPECINDEX_1RE'][hd_indx,44], dapall['SPECINDEX_1RE'][hd_indx,21],
                    lw=0, marker='.', s=10, zorder=2,
                    c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][hd_indx]],
                    norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
                    cmap='viridis')

    cax = fig.add_axes([0.51, 0.48, 0.15, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[1e9,1e10,1e11])
                         #, ticks=[0.1,1], format=ticker.ScalarFormatter())
#    pyplot.setp(cb.ax.xaxis.get_ticklabels(), fontsize=10)
    cax.text(0.5, 2.8, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
             ha='center', va='center', transform=cax.transAxes) #, fontsize=10)

    ax.text(-0.25, 0.5, r'$\langle {\rm H}\delta_{\rm A}\rangle_e$ [$\AA$, in absorption]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax.text(0.5, -0.18, r'$\langle {\rm Dn4000}\rangle_e$',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

#    ax.text(-0.22, 0.5, r'$\Delta V_\ast$ [km/s]',
#            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
#            rotation='vertical')
#    ax.text(0.55, 0.89, r'$\mathcal{M}_\ast\propto\ \Delta V^4$',
#            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
#            rotation=20, zorder=4)

    ofile = '../ms/figs/ew_d4000.pdf'
    fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    dapall_file = '/Users/westfall/Work/MaNGA/testing/dr15/dapall-v2_4_3-2.2.1.fits'
    dapall = fits.open(dapall_file)['DAPALL'].data

    drpall_file = '/Users/westfall/Work/MaNGA/testing/dr15/drpall-v2_4_3.fits'
    drpall = fits.open(drpall_file)[1].data

    ew_d4000(drpall, dapall)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))

