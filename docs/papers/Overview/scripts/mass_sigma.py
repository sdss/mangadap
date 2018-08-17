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
    ax = fig.add_axes(pos, facecolor='0.97')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    return ax


def mass_sigma(dapall, drpall):
    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == 'HYB10-GAU-MILESHC')
    indx &= (dapall['STELLAR_SIGMA_1RE'] > 0)
    indx &= (dapall['NSA_ELPETRO_TH50_R'] > 0)
    indx &= (dapall['ADIST_Z'] > 0)

    rekpc = dapall['NSA_ELPETRO_TH50_R']*numpy.radians(1/3600)*1e3*dapall['ADIST_Z']

    font = { 'size' : 18 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    mass_lim = [7e7, 3e12]
    sigm_lim = [10,  800]

    inst = numpy.array([[mass_lim[0], 70],
                        [mass_lim[1], 70],
                        [mass_lim[1], sigm_lim[0]],
                        [mass_lim[0], sigm_lim[0]]])

    ax = init_ax(fig, [0.2, 0.20, 0.6, 0.6])
    ax.set_xlim(mass_lim)
    ax.set_ylim(sigm_lim)
    ax.set_xscale('log')
    ax.set_yscale('log')

    x = numpy.linspace(6, 14, 100)
#    y = x/4 - 0.25
#    ax.plot(numpy.power(10,x), numpy.power(10,y-1.5), color='0.8', zorder=0)
#    ax.plot(numpy.power(10,x), numpy.power(10,y-0.75), color='0.8', zorder=0)
#    ax.plot(numpy.power(10,x), numpy.power(10,y), color='0.8', zorder=0)
#    ax.plot(numpy.power(10,x), numpy.power(10,y+0.75), color='0.8', zorder=0)

    y = x/2 - 2.6
    ax.plot(numpy.power(10,x), numpy.power(10,y-2), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y-1), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y+1), color='0.8', zorder=0)

    ax.fill(inst[:,0], inst[:,1], color='0.1', alpha=0.1, zorder=1)

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][indx]],
                    dapall['STELLAR_SIGMA_1RE'][indx], lw=0, marker='.', s=20, zorder=3,
                    c=rekpc[indx],
#                    vmin=0.1, vmax=18, cmap='inferno')
                    norm=colors.LogNorm(vmin=1, vmax=30), cmap='inferno')

    cax = fig.add_axes([0.63, 0.26, 0.15, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[1,10],
                         format=ticker.ScalarFormatter())
    cax.text(0.5, 3.0, r'$R_e$ [kpc]', ha='center', va='center', transform=cax.transAxes)
    

    ax.text(0.5, -0.15, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.2, 0.5, r'$\sigma_{\ast,e}$ [km/s]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')
    ax.text(0.51, 0.85, r'$\mathcal{M}_\ast\propto\ \sigma^2_{\ast,e}$',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation=53, zorder=4)
    ax.text(0.99, 0.37, r'$\sigma_{\ast,e}\ \lesssim\ \sigma_{\rm inst}$',
            horizontalalignment='right', verticalalignment='center', transform=ax.transAxes,
            zorder=4)

    ofile = '../ms/figs/mass_sigma.pdf'
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

    mass_sigma(dapall, drpall)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))


