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
    ax = fig.add_axes(pos)#, facecolor='0.98')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    return ax


def velocity_gradient(drpall, dapall):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == 'HYB10-GAU-MILESHC')
    indx &= (drpall['nsa_elpetro_ba'][dapall['DRPALLINDX']] > 0)

    ha_indx = indx & (dapall['HA_GVEL_LO_CLIP'] > -999) & (dapall['HA_GVEL_HI_CLIP'] > -999)
    ha_vdiff = dapall['HA_GVEL_HI_CLIP'][ha_indx] - dapall['HA_GVEL_LO_CLIP'][ha_indx]
    ha_sini = numpy.ma.sqrt(1-numpy.square(drpall['nsa_elpetro_ba'][dapall['DRPALLINDX'][ha_indx]]))

    st_indx = indx & (dapall['STELLAR_VEL_LO_CLIP'] > -999) & (dapall['STELLAR_VEL_HI_CLIP'] > -999)
    st_vdiff = dapall['STELLAR_VEL_HI_CLIP'][st_indx] - dapall['STELLAR_VEL_LO_CLIP'][st_indx]
    st_sini = numpy.ma.sqrt(1-numpy.square(drpall['nsa_elpetro_ba'][dapall['DRPALLINDX'][st_indx]]))

    mass_lim = [2e8, 5e11]
    velw_lim = [30,  4000]

    ha_sb_lim = [0.02, 8]
    st_sb_lim = [0.03, 1]

    font = { 'size' : 14 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.3, 0.56, 0.4, 0.4])
    ax.set_xlim(mass_lim)
    ax.set_ylim(velw_lim)
    ax.tick_params(axis='x', which='both', direction='in')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.NullFormatter())

    x = numpy.linspace(6, 14, 100)
    y = x/4
    ax.plot(numpy.power(10,x), numpy.power(10,y-0.75), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y+0.75), color='0.8', zorder=0)

#    y = 0.25 * x + 0.0
#    ax.plot(numpy.power(10,x), numpy.power(10,y), color='C3', zorder=1)

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][ha_indx]],
               numpy.ma.divide(ha_vdiff, ha_sini),
               c=dapall['EMLINE_GSB_1RE'][ha_indx,18],
               lw=0, marker='.', s=10, zorder=2,
               norm=colors.LogNorm(vmin=ha_sb_lim[0], vmax=ha_sb_lim[1]),
               cmap='inferno_r')
    cax = fig.add_axes([0.54, 0.60, 0.15, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[0.1,1],
                         format=ticker.ScalarFormatter())
    pyplot.setp(cb.ax.xaxis.get_ticklabels(), fontsize=10)
    cax.text(0.8, 2.5, r'$\langle I_{{\rm H}\alpha}\rangle_e$', ha='center', va='center',
             transform=cax.transAxes, fontsize=10)

    ax.text(-0.22, 0.5, r'$\Delta V_{\rm gas}$ [km/s]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')


    ax = init_ax(fig, [0.3, 0.15, 0.4, 0.4])
    ax.set_xlim(mass_lim)
    ax.set_ylim(velw_lim)
    ax.set_xscale('log')
    ax.set_yscale('log')

    x = numpy.linspace(6, 14, 100)
    y = x/4
    ax.plot(numpy.power(10,x), numpy.power(10,y-0.75), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y+0.75), color='0.8', zorder=0)

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][st_indx]],
                    numpy.ma.divide(st_vdiff, st_sini), c=dapall['SB_1RE'][st_indx],
                    lw=0, marker='.', s=10, zorder=2,
                    norm=colors.LogNorm(vmin=st_sb_lim[0], vmax=st_sb_lim[1]),
                    cmap='inferno_r')

    cax = fig.add_axes([0.31, 0.50, 0.11, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[0.1,1],
                         format=ticker.ScalarFormatter())
    pyplot.setp(cb.ax.xaxis.get_ticklabels(), fontsize=10)
    cax.text(0.5, 2.5, r'$\langle I_g\rangle_e$', ha='center', va='center',
             transform=cax.transAxes, fontsize=10)

    ax.text(0.5, -0.18, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.22, 0.5, r'$\Delta V_\ast$ [km/s]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')
    ax.text(0.55, 0.89, r'$\mathcal{M}_\ast\propto\ \Delta V^4$',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation=20, zorder=4)

    ofile = '../ms/figs/mass_vel.pdf'
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

    velocity_gradient(drpall, dapall)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))


