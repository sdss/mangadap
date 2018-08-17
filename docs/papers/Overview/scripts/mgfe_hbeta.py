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

def mgfe_hbeta(drpall, dapall):
    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == 'HYB10-GAU-MILESHC')
    indx &= (drpall['nsa_sersic_mass'][dapall['DRPALLINDX']] > 0)

    indx &= (dapall['SPECINDEX_1RE'][:,12] > -999)
    indx &= (dapall['SPECINDEX_1RE'][:,13] > -999)
    indx &= (dapall['SPECINDEX_1RE'][:,14] > -999)
    indx &= (dapall['SPECINDEX_1RE'][:,8] > -999)

    mgfe = numpy.ma.sqrt(dapall['SPECINDEX_1RE'][indx,12] * (0.72*dapall['SPECINDEX_1RE'][indx,13] 
                        + 0.28*dapall['SPECINDEX_1RE'][indx,13]))

    mass_lim = [1e9, 7e11]
    mgfe_lim = [0, 4.5]
    hb_lim = [0.5, 5.9]

    font = { 'size' : 16 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.2, 0.2, 0.6, 0.6])
    ax.set_xlim(mgfe_lim)
    ax.set_ylim(hb_lim)
    ax.tick_params(axis='x', which='both', direction='in')

    sc = ax.scatter(mgfe, dapall['SPECINDEX_1RE'][indx,8],
                    marker='.', s=30, lw=0, zorder=2,
                    c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][indx]],
                    norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
                    cmap='viridis')

    ax.text(-0.15, 0.5, r'$\langle {\rm H}\beta\rangle_e$ [$\AA$, in absorption]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax.text(0.5, -0.12, r'[$\langle$Mgb$\rangle_e\cdot$(0.72 $\langle$Fe5270$\rangle_e$ + 0.28 $\langle$Fe5335$\rangle_e$)]$^{1/2}$',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.5, -0.19, r'[$\AA$, in absorption]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    cax = fig.add_axes([0.58, 0.71, 0.20, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[1e9,1e10,1e11])
    cax.text(0.5, 3.0, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
             ha='center', va='center', transform=cax.transAxes) #, fontsize=10)

    ofile = '../ms/figs/mgfe_hbeta.pdf'
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

    mgfe_hbeta(drpall, dapall)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))
