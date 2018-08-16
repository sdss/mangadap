#!/usr/bin/env python3

import os
import time
import numpy

from matplotlib import pyplot, ticker, rc

from astropy.io import fits
import astropy.constants

from mangadap.proc.util import growth_lim
from mangadap.util.bitmask import BitMask

#-----------------------------------------------------------------------------
def init_ax(fig, pos):
    ax = fig.add_axes(pos, facecolor='0.98')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    return ax


def redshift_distribution(dapall):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == 'HYB10-GAU-MILESHC')
    nsa_indx = indx & (dapall['NSA_Z'] > -999)
    str_indx = indx & (dapall['STELLAR_Z'] > -999)
    gas_indx = indx & (dapall['HA_Z'] > -999)

    font = { 'size' : 14 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    offset_lim = [-90, 90]
    ax = init_ax(fig, [0.17, 0.10, 0.56, 0.3])
    ax.set_xlim([-0.01, 0.16])
    ax.set_ylim(offset_lim)
    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')

    dcz = astropy.constants.c.to('km/s').value * (dapall['STELLAR_Z'][str_indx & nsa_indx]
                                                  - dapall['NSA_Z'][str_indx & nsa_indx])
    ax.scatter(dapall['NSA_Z'][str_indx & nsa_indx], dcz,
                   marker='.', s=10, lw=0, color='k', zorder=3)
    ax.text(0.5, -0.2, r'$z_{\rm NSA}$', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.18, 0.5, r'$c (z_\ast - z_{\rm NSA})$ [km/s]', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, rotation='vertical')

    ax = init_ax(fig, [0.73, 0.10, 0.1, 0.3])
    ax.set_ylim(offset_lim)
    ax.tick_params(which='both', direction='in', bottom=False)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.grid(True, which='major', axis='y', color='0.8', zorder=0, linestyle='-')
    ax.hist(dcz, range=offset_lim, bins=100, alpha=0.7, color='0.4', histtype='stepfilled',
            orientation='horizontal', zorder=3)
    ax.hist(dcz, range=offset_lim, bins=100, alpha=0.7, color='k', histtype='step', lw=0.5,
            orientation='horizontal', zorder=4)

    ax = init_ax(fig, [0.22, 0.50, 0.56, 0.3])
    ax.set_xlim(offset_lim)
    ax.tick_params(which='both', left=False)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.grid(True, which='major', axis='x', color='0.8', zorder=0, linestyle='-')
    dcz = astropy.constants.c.to('km/s').value * (dapall['STELLAR_Z'][gas_indx & str_indx]
                                                  - dapall['HA_Z'][gas_indx & str_indx]),
    ax.hist(dcz, range=offset_lim, bins=100, alpha=0.7, color='0.4', histtype='stepfilled',
            zorder=3)
    ax.hist(dcz, range=offset_lim, bins=100, alpha=0.7, color='k', histtype='step', lw=0.5,
            zorder=4)
    ax.text(0.5, -0.2, r'$c (z_\ast - z_{\rm gas})$ [km/s]', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes)

    ofile = '../ms/figs/redshift_difference.pdf'
    fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    dapall_file = '/Users/westfall/Work/MaNGA/testing/dr15/dapall-v2_4_3-2.2.1.fits'
    dapall = fits.open(dapall_file)['DAPALL'].data
    redshift_distribution(dapall)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))


