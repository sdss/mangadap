import os
import time
import argparse
import warnings

from IPython import embed

import numpy

from matplotlib import pyplot, ticker, rc, colors

from astropy.io import fits
import astropy.constants

from mangadap.util.bitmask import BitMask
from mangadap.dapfits import DAPQualityBitMask
from mangadap.config import defaults
from mangadap.proc.util import growth_lim
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.util.fileio import channel_dictionary

from mangadap.scripts import scriptbase

#-----------------------------------------------------------------------------
def init_ax(fig, pos, facecolor=None, grid=False):
    ax = fig.add_axes(pos, facecolor=facecolor)
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    if grid:
        ax.grid(True, which='major', color='0.8', zorder=1, linestyle='-', lw=0.5)
    return ax


def dapall_qa(drpall, dapall, analysis_path, daptype, eml, spi):
    oroot = os.path.join(analysis_path, daptype, 'qa')
    if not os.path.isdir(oroot):
        os.makedirs(oroot)

    # Plot H-delta A equivalent width versus D4000
    ew_d4000(drpall, dapall, daptype, eml, spi, ofile=os.path.join(oroot, 'dapall_ew_d4000.png'))

    # Plot the star-formation, mass diagram
    mass_lha(drpall, dapall, daptype, eml, ofile=os.path.join(oroot, 'dapall_mass_lha.png'))

    # Plot the approximate Faber-Jackson relation
    mass_sigma(drpall, dapall, daptype, ofile=os.path.join(oroot, 'dapall_mass_sigma.png'))

    # Plot the approximate Tully-Fisher relation
    velocity_gradient(drpall, dapall, daptype, eml,
                      ofile=os.path.join(oroot, 'dapall_mass_vel.png'))

    # Plot the MgFe and Hbeta absorption indices
    mgfe_hbeta(drpall, dapall, daptype, spi, ofile=os.path.join(oroot, 'dapall_mgfe_hbeta.png'))

    # Plot the radial coverage of each cube
    radial_coverage_histogram(dapall, daptype,
                              ofile=os.path.join(oroot, 'dapall_radialcoverage.png'))

    # Plot the redshift differences
    redshift_distribution(dapall, daptype, ofile=os.path.join(oroot, 'dapall_redshift_dist.png'))


def ew_d4000(drpall, dapall, daptype, eml, spi, ofile=None):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype)
    indx &= (drpall['nsa_sersic_mass'][dapall['DRPALLINDX']] > 0)

    ha_indx = indx & (dapall['EMLINE_GEW_1RE'][:,eml['Ha-6564']] > -999) \
                    & (dapall['SPECINDEX_1RE'][:,spi['Dn4000']] > -999)
    hd_indx = indx & (dapall['SPECINDEX_1RE'][:,spi['HDeltaA']] > -999) \
                    & (dapall['SPECINDEX_1RE'][:,spi['Dn4000']] > -999)

    is_critical = DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL')

    ha_bd = ha_indx & is_critical
    hd_bd = hd_indx & is_critical

    ha_gd = ha_indx & numpy.invert(is_critical)
    hd_gd = hd_indx & numpy.invert(is_critical)

    mass_lim = [1e9, 7e11]
    d4_lim = [0.9, 2.3]
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

    sc = ax.scatter(dapall['SPECINDEX_1RE'][ha_gd,spi['Dn4000']],
                    dapall['EMLINE_GEW_1RE'][ha_gd,eml['Ha-6564']],
                    lw=0, marker='.', s=20, zorder=3,
                    c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][ha_gd]],
                    norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
                    cmap='viridis')

    ax.scatter(dapall['SPECINDEX_1RE'][ha_bd,spi['Dn4000']],
               dapall['EMLINE_GEW_1RE'][ha_bd,eml['Ha-6564']],
               lw=0, marker='x', s=10, zorder=2,
               c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][ha_bd]],
               norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
               cmap='viridis')

    ax.text(-0.25, 0.5, r'$\langle {\rm EW}_{{\rm H}\alpha}\rangle_e$ [$\AA$, in emission]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax = init_ax(fig, [0.3, 0.15, 0.4, 0.4])
    ax.set_xlim(d4_lim)
    ax.set_ylim(hd_lim)

    sc = ax.scatter(dapall['SPECINDEX_1RE'][hd_gd,spi['Dn4000']],
                    dapall['SPECINDEX_1RE'][hd_gd,spi['HDeltaA']],
                    lw=0, marker='.', s=20, zorder=3,
                    c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][hd_gd]],
                    norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
                    cmap='viridis')

    ax.scatter(dapall['SPECINDEX_1RE'][hd_bd,spi['Dn4000']],
               dapall['SPECINDEX_1RE'][hd_bd,spi['HDeltaA']],
               lw=0, marker='x', s=10, zorder=2,
               c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][hd_bd]],
               norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
               cmap='viridis')

    cax = fig.add_axes([0.51, 0.48, 0.15, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[1e9,1e10,1e11])
                         #, ticks=[0.1,1], format=ticker.ScalarFormatter())
    cax.text(0.5, 2.8, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
             ha='center', va='center', transform=cax.transAxes) #, fontsize=10)

    ax.text(-0.25, 0.5, r'$\langle {\rm H}\delta_{\rm A}\rangle_e$ [$\AA$, in absorption]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax.text(0.5, -0.18, r'$\langle {\rm Dn4000}\rangle_e$',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def mass_lha(drpall, dapall, daptype, eml, ofile=None):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype)
    indx &= (drpall['nsa_sersic_mass'][dapall['DRPALLINDX']] > 0)

    indx &= (dapall['EMLINE_GSB_1RE'][:,eml['Hb-4862']] > -999)
    indx &= (dapall['EMLINE_GSB_1RE'][:,eml['Ha-6564']] > -999)
    indx &= (dapall['SFR_1RE'] > -999)

    decrement = numpy.ma.divide(dapall['EMLINE_GSB_1RE'][:,eml['Ha-6564']],
                                   dapall['EMLINE_GSB_1RE'][:,eml['Hb-4862']])

    loew_indx =  indx & (dapall['EMLINE_GEW_1RE'][:,eml['Ha-6564']] > -999)\
                        & (dapall['EMLINE_GEW_1RE'][:,eml['Ha-6564']] < 2)
    hiew_indx = indx & (dapall['EMLINE_GEW_1RE'][:,eml['Ha-6564']] > 2)

    is_critical = DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL')

    loew_bd =  loew_indx & is_critical
    hiew_bd =  hiew_indx & is_critical

    loew_gd =  loew_indx & numpy.invert(is_critical)
    hiew_gd =  hiew_indx & numpy.invert(is_critical)


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

    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][loew_gd]],
               numpy.ma.log10(dapall['SFR_1RE'][loew_gd])+41.27,
               color='0.5', marker='.', s=20, lw=0, zorder=3, alpha=0.3)

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][hiew_gd]],
                    numpy.ma.log10(dapall['SFR_1RE'][hiew_gd])+41.27,
                    c=decrement[hiew_gd], vmin=2.8, vmax=7,
                    marker='.', s=20, lw=0, zorder=4, cmap='inferno')

    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][loew_bd]],
               numpy.ma.log10(dapall['SFR_1RE'][loew_bd])+41.27,
               color='0.5', marker='x', s=10, lw=0, zorder=3, alpha=0.3)

    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][hiew_bd]],
               numpy.ma.log10(dapall['SFR_1RE'][hiew_bd])+41.27,
               c=decrement[hiew_bd], vmin=2.8, vmax=7,
               marker='x', s=10, lw=0, zorder=4, cmap='inferno')

    ax.text(-0.15, 0.5, r'$\log\ \langle L_{{\rm H}\alpha}\rangle_e$ [erg/s]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
            rotation='vertical')

    ax.text(0.5, -0.12, r'$\mathcal{M}_\ast$ [$h^{-2} \mathcal{M}_\odot$]',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    cax = fig.add_axes([0.58, 0.25, 0.20, 0.01])
    cb = pyplot.colorbar(sc, cax=cax, orientation='horizontal', ticks=[3, 4, 5, 6, 7])
    cax.text(0.5, 3.0, r'$L_{{\rm H}\alpha}/L_{{\rm H}\beta}$',
             ha='center', va='center', transform=cax.transAxes) #, fontsize=10)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def mass_sigma(drpall, dapall, daptype, ofile=None):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype)
    indx &= (dapall['STELLAR_SIGMA_1RE'] > 0)
    indx &= (dapall['NSA_ELPETRO_TH50_R'] > 0)
    indx &= (dapall['ADIST_Z'] > 0)

    is_critical = DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL')
    bd =  indx & is_critical
    gd =  indx & numpy.invert(is_critical)

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
    y = x/2 - 2.6
    ax.plot(numpy.power(10,x), numpy.power(10,y-2), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y-1), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y), color='0.8', zorder=0)
    ax.plot(numpy.power(10,x), numpy.power(10,y+1), color='0.8', zorder=0)

    ax.fill(inst[:,0], inst[:,1], color='0.1', alpha=0.1, zorder=1)

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][gd]],
                    dapall['STELLAR_SIGMA_1RE'][gd], lw=0, marker='.', s=20, zorder=3,
                    c=rekpc[gd], norm=colors.LogNorm(vmin=1, vmax=30), cmap='inferno')

    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][bd]],
               dapall['STELLAR_SIGMA_1RE'][bd], lw=0, marker='x', s=10, zorder=3,
               c=rekpc[bd], norm=colors.LogNorm(vmin=1, vmax=30), cmap='inferno')

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

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def velocity_gradient(drpall, dapall, daptype, eml, ofile=None):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype)
    indx &= (drpall['nsa_elpetro_ba'][dapall['DRPALLINDX']] > 0)

    ha_indx = indx & (dapall['HA_GVEL_LO_CLIP'] > -999) & (dapall['HA_GVEL_HI_CLIP'] > -999)

    ha_vdiff = dapall['HA_GVEL_HI_CLIP'] - dapall['HA_GVEL_LO_CLIP']
    ha_sini = numpy.ma.sqrt(1-numpy.square(drpall['nsa_elpetro_ba'][dapall['DRPALLINDX']]))

    st_indx = indx & (dapall['STELLAR_VEL_LO_CLIP'] > -999) & (dapall['STELLAR_VEL_HI_CLIP'] > -999)

    st_vdiff = dapall['STELLAR_VEL_HI_CLIP'] - dapall['STELLAR_VEL_LO_CLIP']
    st_sini = numpy.ma.sqrt(1-numpy.square(drpall['nsa_elpetro_ba'][dapall['DRPALLINDX']]))

    is_critical = DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL')

    ha_bd =  ha_indx & is_critical
    st_bd =  st_indx & is_critical

    ha_gd =  ha_indx & numpy.invert(is_critical)
    st_gd =  st_indx & numpy.invert(is_critical)

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

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][ha_gd]],
               numpy.ma.divide(ha_vdiff[ha_gd], ha_sini[ha_gd]),
               c=dapall['EMLINE_GSB_1RE'][ha_gd,eml['Ha-6564']],
               lw=0, marker='.', s=20, zorder=2,
               norm=colors.LogNorm(vmin=ha_sb_lim[0], vmax=ha_sb_lim[1]),
               cmap='inferno_r')
    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][ha_bd]],
               numpy.ma.divide(ha_vdiff[ha_bd], ha_sini[ha_bd]),
               c=dapall['EMLINE_GSB_1RE'][ha_bd,eml['Ha-6564']],
               lw=0, marker='x', s=10, zorder=2,
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

    sc = ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][st_gd]],
                    numpy.ma.divide(st_vdiff[st_gd], st_sini[st_gd]), c=dapall['SB_1RE'][st_gd],
                    lw=0, marker='.', s=20, zorder=2,
                    norm=colors.LogNorm(vmin=st_sb_lim[0], vmax=st_sb_lim[1]),
                    cmap='inferno_r')
    ax.scatter(drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][st_bd]],
               numpy.ma.divide(st_vdiff[st_bd], st_sini[st_bd]), c=dapall['SB_1RE'][st_bd],
               lw=0, marker='x', s=10, zorder=2,
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

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def mgfe_hbeta(drpall, dapall, daptype, spi, ofile=None):
    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype)
    indx &= (drpall['nsa_sersic_mass'][dapall['DRPALLINDX']] > 0)

    indx &= (dapall['SPECINDEX_1RE'][:,spi['Mgb']] > -999)
    indx &= (dapall['SPECINDEX_1RE'][:,spi['Fe5270']] > -999)
    indx &= (dapall['SPECINDEX_1RE'][:,spi['Fe5335']] > -999)
    indx &= (dapall['SPECINDEX_1RE'][:,spi['Hb']] > -999)

    is_critical = DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL')
    bd =  indx & is_critical
    gd =  indx & numpy.invert(is_critical)

    mgfe = numpy.ma.sqrt(dapall['SPECINDEX_1RE'][:,spi['Mgb']]
                            * (0.72*dapall['SPECINDEX_1RE'][:,spi['Fe5270']] 
                                + 0.28*dapall['SPECINDEX_1RE'][:,spi['Fe5335']]))

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

    sc = ax.scatter(mgfe[gd], dapall['SPECINDEX_1RE'][gd,spi['Hb']],
                    marker='.', s=30, lw=0, zorder=2,
                    c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][gd]],
                    norm=colors.LogNorm(vmin=mass_lim[0], vmax=mass_lim[1]),
                    cmap='viridis')
    ax.scatter(mgfe[bd], dapall['SPECINDEX_1RE'][bd,spi['Hb']],
               marker='x', s=15, lw=0, zorder=2,
               c=drpall['nsa_sersic_mass'][dapall['DRPALLINDX'][bd]],
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

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


def radial_coverage_histogram(dapall, daptype, ofile=None):

    # Must have finished, be of the correct daptype, and *not* be
    # critical
    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype) \
                & numpy.invert(DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL'))

    # Get the range
    arcsec_rng = growth_lim(dapall['RCOV90'][indx], 0.99, fac=1.1)
    rore = numpy.ma.divide(dapall['RCOV90'][indx], dapall['NSA_ELPETRO_TH50_R'][indx])
    rore[dapall['NSA_ELPETRO_TH50_R'][indx] < 0] = numpy.ma.masked
    rore_rng = growth_lim(rore, 0.99, fac=1.1)

    # Determine the sample
    sdssMaskbits = defaults.sdss_maskbits_file()
    targ1bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET1')

    ancillary = dapall['MNGTARG3'][indx] > 0
    primaryplus = targ1bm.flagged(dapall['MNGTARG1'][indx],
                                  flag=['PRIMARY_v1_2_0', 'COLOR_ENHANCED_v1_2_0'])
    secondary = targ1bm.flagged(dapall['MNGTARG1'][indx], flag='SECONDARY_v1_2_0')
    other = numpy.invert(ancillary) & numpy.invert(primaryplus) & numpy.invert(secondary)

    print('Num. of observations with 90% coverage <1 arcsec: {0}'.format(
            numpy.sum(dapall['RCOV90'][indx] < 1)))

    font = { 'size' : 10 }
    rc('font', **font)

    # Instantiate the figure
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    nlim = [0.8, 1000]

    # Show RCOV90 in terms of arcseconds for each IFU
    ax = init_ax(fig, [0.1, 0.58, 0.8, 0.4])
    ax.set_ylim(nlim)
    ax.set_yscale('log', nonpositive='clip')
    ifu = [ 19, 37, 61, 91, 127 ]
    for i,f in enumerate(ifu):
        _indx = indx & ((dapall['IFUDESIGN']/100).astype(int) == f)
        if not numpy.any(_indx):
            continue
        ax.hist(dapall['RCOV90'][_indx], bins=50, range=arcsec_rng, alpha=0.5,
                histtype='stepfilled', color='C{0}'.format(i), zorder=2)
        ax.hist(dapall['RCOV90'][_indx], bins=50, range=arcsec_rng,
                histtype='step', color='C{0}'.format(i), zorder=3)
        medp = numpy.median(dapall['RCOV90'][_indx])
        ax.plot([medp, medp], nlim, color='C{0}'.format(i), zorder=4, linestyle='--')
    ax.text(-0.1, 0.5, 'N', ha='center', va='center', transform=ax.transAxes, rotation='vertical')
    ax.text(0.5, -0.13, r'$R_{90}$ (arcsec)', ha='center', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.95, '19-fiber', color='C0', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.90, '37-fiber', color='C1', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.85, '61-fiber', color='C2', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.80, '91-fiber', color='C3', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.75, '127-fiber', color='C4', ha='left', va='center', transform=ax.transAxes)

    # Show RCOV90 in terms of Re for each sample
    ax = init_ax(fig, [0.1, 0.08, 0.8, 0.4])
    ax.set_ylim(nlim)
    ax.set_yscale('log', nonpositive='clip')

    if numpy.any(primaryplus):
        _rore = (rore[primaryplus]).compressed()
        ax.hist(_rore, bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
                color='C0', zorder=2)
        ax.hist(_rore, bins=50, range=rore_rng, histtype='step',
                color='C0', zorder=3)
        medp = numpy.median(_rore)
        ax.plot([medp, medp], nlim, color='C0', zorder=4, linestyle='--')

    if numpy.any(secondary):
        _rore = (rore[secondary]).compressed()
        ax.hist(_rore, bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
                color='C1', zorder=2)
        ax.hist(_rore, bins=50, range=rore_rng, histtype='step',
                color='C1', zorder=3)
        medp = numpy.median(_rore)
        ax.plot([medp, medp], nlim, color='C1', zorder=4, linestyle='--')

    if numpy.any(other):
        _rore = (rore[other]).compressed()
        ax.hist(_rore, bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
                color='C2', zorder=2)
        ax.hist(_rore, bins=50, range=rore_rng, histtype='step',
                color='C2', zorder=3)

    if numpy.any(ancillary):
        _rore = (rore[ancillary]).compressed()
        ax.hist(_rore, bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
                color='C3', zorder=2)
        ax.hist(_rore, bins=50, range=rore_rng, histtype='step',
                color='C3', zorder=3)

    ax.text(-0.1, 0.5, 'N', ha='center', va='center', transform=ax.transAxes, rotation='vertical')
    ax.text(0.5, -0.12, r'$R_{90}/R_{eff}$', ha='center', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.95, 'Primary+', color='C0', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.90, 'Secondary', color='C1', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.85, 'Other', color='C2', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.80, 'Ancillary', color='C3', ha='left', va='center', transform=ax.transAxes)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


def redshift_distribution(dapall, daptype, ofile=None):

    indx = (dapall['DAPDONE'] == 1) & (dapall['DAPTYPE'] == daptype)
    ns_indx = indx & (dapall['NSA_Z'] > -999) & (dapall['STELLAR_Z'] > -999)
    gs_indx = indx & (dapall['HA_Z'] > -999) & (dapall['STELLAR_Z'] > -999)

    is_critical = DAPQualityBitMask().flagged(dapall['DAPQUAL'], 'CRITICAL')

    ns_bd = ns_indx & is_critical
    gs_bd = gs_indx & is_critical

    ns_gd = ns_indx & numpy.invert(is_critical)
    gs_gd = gs_indx & numpy.invert(is_critical)


    font = { 'size' : 14 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    offset_lim = [-90, 90]
    ax = init_ax(fig, [0.17, 0.10, 0.56, 0.3])
    ax.set_xlim([-0.01, 0.16])
    ax.set_ylim(offset_lim)
    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')

    dcz = astropy.constants.c.to('km/s').value * (dapall['STELLAR_Z'] - dapall['NSA_Z'])

    ax.scatter(dapall['NSA_Z'][ns_gd], dcz[ns_gd], marker='.', s=20, lw=0, color='k', zorder=3)
    ax.scatter(dapall['NSA_Z'][ns_bd], dcz[ns_bd], marker='.', s=20, lw=0, color='C3', zorder=3)
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
    ax.hist(dcz[ns_gd], range=offset_lim, bins=100, alpha=0.7, color='0.4', histtype='stepfilled',
            orientation='horizontal', zorder=3)
    ax.hist(dcz[ns_gd], range=offset_lim, bins=100, alpha=0.7, color='k', histtype='step', lw=0.5,
            orientation='horizontal', zorder=4)

    ax = init_ax(fig, [0.22, 0.50, 0.56, 0.3])
    ax.set_xlim(offset_lim)
    ax.tick_params(which='both', left=False)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.grid(True, which='major', axis='x', color='0.8', zorder=0, linestyle='-')

    dcz = astropy.constants.c.to('km/s').value * (dapall['STELLAR_Z'] - dapall['HA_Z'])

    ax.hist(dcz[gs_gd], range=offset_lim, bins=100, alpha=0.7, color='0.4', histtype='stepfilled',
            zorder=3)
    ax.hist(dcz[gs_gd], range=offset_lim, bins=100, alpha=0.7, color='k', histtype='step', lw=0.5,
            zorder=4)
    ax.text(0.5, -0.2, r'$c (z_\ast - z_{\rm gas})$ [km/s]', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


class DapAllQA(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'dapall_qa'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Construct QA plots using the DAPall metadata',
                                    width=width)

        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument("--redux_path", type=str, help="main DRP output path", default=None)
        parser.add_argument("--analysis_path", type=str, help="main DAP output path", default=None)

        parser.add_argument("--plan_file", type=str, help="parameter file with the MaNGA DAP "
                            "execution plan to use instead of the default" , default=None)

        parser.add_argument('--daptype', type=str, help='DAP processing type', default=None)
        parser.add_argument('--normal_backend', dest='bgagg', action='store_false', default=True)

        return parser

    @staticmethod
    def main(args):
        t = time.perf_counter()

        if args.bgagg:
            pyplot.switch_backend('agg')

        # Set the paths
        redux_path = defaults.drp_redux_path(drpver=args.drpver) \
                            if args.redux_path is None else args.redux_path
        analysis_path = defaults.dap_analysis_path(drpver=args.drpver, dapver=args.dapver) \
                                if args.analysis_path is None else args.analysis_path

        daptypes = []
        if args.daptype is None:
            analysisplan = AnalysisPlanSet.default() if args.plan_file is None \
                                else AnalysisPlanSet.from_par_file(args.plan_file)
            for p in analysisplan:
                bin_method = SpatiallyBinnedSpectra.define_method(p['bin_key'])
                sc_method = StellarContinuumModel.define_method(p['continuum_key'])
                el_method = EmissionLineModel.define_method(p['elfit_key'])
                daptypes += [defaults.dap_method(bin_method['key'],
                                                 sc_method['fitpar']['template_library_key'],
                                                 el_method['continuum_tpl_key'])]
        else:
            daptypes = [args.daptype]

        drpall_file = defaults.drpall_file(drpver=args.drpver, redux_path=redux_path)
        dapall_file = defaults.dapall_file(drpver=args.drpver, dapver=args.dapver,
                                           analysis_path=analysis_path)

        drpall_hdu = fits.open(drpall_file)
        dapall_hdu = fits.open(dapall_file)

        eml = channel_dictionary(dapall_hdu, 0, prefix='ELG')
        spi = channel_dictionary(dapall_hdu, 0, prefix='SPI')

        for daptype in daptypes:
            dapall_qa(drpall_hdu['MANGA'].data, dapall_hdu[daptype].data, analysis_path,
                      daptype, eml, spi)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

