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
    ax = fig.add_axes(pos, facecolor='0.98')
    ax.grid(True, which='major', color='0.9', zorder=0, linestyle='-')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    return ax

def mass_lha(drpall, dapall):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == 'HYB10-GAU-MILESHC')
    indx &= (drpall['nsa_sersic_mass'][dapall['DRPALLINDX']] > 0)

    indx &= (dapall['EMLINE_GSB_1RE'][:,11] > -999)
    indx &= (dapall['EMLINE_GSB_1RE'][:,18] > -999)
    indx &= (dapall['SFR_1RE'] > -999)

    decrement = numpy.ma.divide(dapall['EMLINE_GSB_1RE'][:,18],
                                   dapall['EMLINE_GSB_1RE'][:,11])

    loew_indx =  indx & (dapall['EMLINE_GEW_1RE'][:,18] > -999)\
                        & (dapall['EMLINE_GEW_1RE'][:,18] < 2)
    hiew_indx = indx & (dapall['EMLINE_GEW_1RE'][:,18] > 2)

    mass_lim = [1e8, 1e12]
    lha_lim = [36.5, 42.5]

    font = { 'size' : 16 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.2, 0.2, 0.6, 0.6])
    ax.set_xlim(mass_lim)
    ax.set_ylim(lha_lim)
    ax.set_xscale('log')

    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][loew_indx]],
               numpy.ma.log10(dapall['SFR_1RE'][loew_indx])+41.27,
               color='0.5', marker='.', s=20, lw=0, zorder=3, alpha=0.3)

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][hiew_indx]],
                    numpy.ma.log10(dapall['SFR_1RE'][hiew_indx])+41.27,
                    c=decrement[hiew_indx], vmin=2.8, vmax=7,
                    marker='.', s=20, lw=0, zorder=4, cmap='inferno')

    ax.text(-0.15, 0.5, r'$\log\ \langle L_{{\rm H}\alpha}\rangle_e$ [erg/s]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax.text(0.5, -0.12, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    cax = fig.add_axes([0.58, 0.25, 0.20, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[3, 4, 5, 6, 7])
    cax.text(0.5, 3.0, r'$L_{{\rm H}\alpha}/L_{{\rm H}\beta}$',
             ha='center', va='center', transform=cax.transAxes) #, fontsize=10)

    ofile = '../ms/figs/mass_lha.pdf'
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

    mass_lha(drpall, dapall)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))

