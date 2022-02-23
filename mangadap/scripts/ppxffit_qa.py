import os
import time

from IPython import embed

import numpy

import argparse

from matplotlib import pyplot, colors, rc, colorbar, ticker, cm

from astropy.io import fits

from mangadap.datacube import MaNGADataCube
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.reductionassessments import ReductionAssessment
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel
from mangadap.proc.util import growth_lim
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.mapping import map_extent, map_beam_patch
from mangadap.config import defaults

from mangadap.scripts import scriptbase


def init_image_ax(fig, pos):
    ax = fig.add_axes(pos, facecolor='0.9')
    ax.minorticks_on()
    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6, direction='in')
    ax.tick_params(which='minor', length=3, direction='in')
    return ax


# THIS IS A TEST FUNCTION!
def stellar_continuum_scatter(plt, ifu, daptype, binid, snr, r68, r99, rchi2,
                              signal, a, da, an, dan, gmr, t, dt, tn, dtn,
                              svel, ssigo, ssigcor, ssigc, extent=None, ofile=None):

    snr_lim = numpy.power(10., growth_lim(numpy.ma.log10(snr), 1.00, fac=1.05))
    snr_lim[0] = max(0.1, snr_lim[0])
    r68_lim = numpy.power(10., growth_lim(numpy.ma.log10(r68), 0.90, fac=1.05))
    r99_lim = numpy.power(10., growth_lim(numpy.ma.log10(r99), 0.90, fac=1.05))
    chi_lim = growth_lim(rchi2, 0.90, fac=1.05)

    s_lim = numpy.power(10., growth_lim(numpy.ma.log10(signal), 1.00, fac=1.10))
    a_lim = growth_lim(numpy.ma.log10(a), 0.99, fac=1.10)
    da_lim = growth_lim(numpy.ma.log10(da), 0.99, fac=1.10)
    an_lim = growth_lim(numpy.ma.log10(an), 1.00, fac=1.10)
    dan_lim = growth_lim(numpy.ma.log10(dan), 1.00, fac=1.10)
    an_lim = numpy.power(10, an_lim)
    dan_lim = numpy.power(10, dan_lim)

    gmr_lim = None if gmr is None else growth_lim(gmr, 0.95, fac=1.10)
    t_lim = growth_lim(numpy.ma.log10(t), 0.99, fac=1.10)
    dt_lim = growth_lim(numpy.ma.log10(dt), 0.99, fac=1.10)
    tn_lim = growth_lim(numpy.ma.log10(tn), 1.00, fac=1.10)
    dtn_lim = growth_lim(numpy.ma.log10(dtn), 1.00, fac=1.10)
    dtn_lim = numpy.power(10, dtn_lim)
    tn_lim = numpy.power(10, tn_lim)

    sv_lim = growth_lim(svel, 0.90, fac=1.2, midpoint=0)
    so_lim = numpy.power(10., growth_lim(numpy.ma.log10(ssigo), 1.00, fac=1.05))
    sc_lim = growth_lim(ssigcor, 0.90, fac=1.05)
    ss_lim = numpy.power(10., growth_lim(numpy.ma.log10(ssigc), 0.90, fac=1.05))

#    so_lim = [ min(so_lim[0],ss_lim[0]), max(so_lim[1],ss_lim[1]) ]
#    ss_lim = so_lim.copy()


    uniq, indx = numpy.unique(binid, return_index=True)
    if uniq[0] == -1:
        indx = indx[1:]

#    fig = pyplot.figure()
#    ax = fig.add_subplot(111, projection='3d')
#
#    ax.scatter(numpy.ma.log10(dtn[indx]), numpy.ma.log10(dan[indx]), numpy.ma.log10(ssigo[indx]),
#               c=snr[indx], cmap='viridis', norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]))
#
#    ax.set_xlim(numpy.log10(dtn_lim))
#    ax.set_ylim(numpy.log10(dan_lim))
#    ax.set_zlim(numpy.log10(so_lim))


#    sc = pyplot.scatter(dan[indx], ssigo[indx], marker='.', s=20, c=rchi2[indx],
#                        cmap='viridis', vmin=chi_lim[0], vmax=chi_lim[1])
#    sc = pyplot.scatter(snr[indx], ssigo[indx], marker='.', s=20, c=dan[indx],
#                        cmap='viridis', norm=colors.LogNorm(vmin=numpy.amin(dan[indx]), vmax=numpy.amax(dan[indx])))
#    sc = pyplot.scatter(signal[indx], dan[indx], marker='.', s=20, c=snr[indx],
#                        cmap='viridis', norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]))
    sc = pyplot.scatter(snr[indx], dan[indx], marker='.', s=20, c=ssigo[indx],
                        cmap='viridis', norm=colors.LogNorm(vmin=so_lim[0], vmax=so_lim[1]))
#    sc = pyplot.scatter(dan[indx], ssigo[indx], marker='.', s=20, c=snr[indx],
#                        cmap='viridis', norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]))
#    sc = pyplot.scatter(dtn[indx], dan[indx], marker='.', s=20, c=snr[indx],
#                        cmap='viridis', norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]))
    cb = pyplot.colorbar(sc)
##    pyplot.xlim(s_lim)
    pyplot.ylim(dan_lim)
    pyplot.xlim(snr_lim)
    pyplot.xscale('log')
    pyplot.yscale('log')

    pyplot.show()
    exit()
    

def masked_imshow(fig, ax, cax, data, extent=None, norm=None, vmin=None, vmax=None, cmap='viridis',
                  zorder=0, cbformat=None, subs=None):
    if numpy.sum(data.mask) != numpy.prod(data.shape):
        if norm is None:
            img = ax.imshow(data, origin='lower', interpolation='nearest', extent=extent,
                            vmin=vmin, vmax=vmax, cmap=cmap, zorder=zorder)
        else:
            img = ax.imshow(data, origin='lower', interpolation='nearest', extent=extent,
                            norm=norm, cmap=cmap, zorder=zorder)

        cb = fig.colorbar(img, cax=cax) if cbformat is None else \
                    fig.colorbar(img, cax=cax, format=cbformat)
        if subs is not None:
            cb.locator = ticker.LogLocator(base=10, subs=(1.,2.,4.,))
            cb.update_ticks()
    else:
        _norm = colors.Normalize(vmin=vmin, vmax=vmax) if norm is None else norm
        cb = colorbar.ColorbarBase(cax, cmap=cm.get_cmap(cmap), norm=_norm)


def stellar_continuum_maps(plt, ifu, daptype, snr, r68, r99, rchi2, signal, a, da, an, dan,
                           gmr, t, dt, tn, dtn, svel, ssigo, ssigcor, ssigc, extent=None,
                           ofile=None):

    font = { 'size' : 6 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(2*w,2*h))

#    Dx = (extent[0]-extent[1])*1.1
#    x_lim = (extent[0]+extent[1])/2 + Dx*numpy.array([1,-1])/2
#    Dy = (extent[2]-extent[3])*1.1
#    y_lim = (extent[2]+extent[3])/2 + Dy*numpy.array([-1,1])/2
#    ax.set_xlim(x_lim)
#    ax.set_ylim(y_lim)

    snr_lim = numpy.power(10., growth_lim(numpy.ma.log10(snr), 0.90, fac=1.05))
    snr_lim[0] = max(0.1, snr_lim[0])
    r68_lim = numpy.power(10., growth_lim(numpy.ma.log10(r68), 0.90, fac=1.05))
    r99_lim = numpy.power(10., growth_lim(numpy.ma.log10(r99), 0.90, fac=1.05))
    chi_lim = growth_lim(rchi2, 0.90, fac=1.05)

    s_lim = numpy.power(10., growth_lim(numpy.ma.log10(signal), 0.99, fac=1.05))
    a_lim = growth_lim(numpy.ma.log10(a), 0.99, fac=1.05)
    da_lim = growth_lim(numpy.ma.log10(da), 0.99, fac=1.05)
    an_lim = growth_lim(numpy.ma.log10(an), 0.99, fac=1.05)
    dan_lim = growth_lim(numpy.ma.log10(dan), 0.99, fac=1.05)

    gmr_lim = None if gmr is None else growth_lim(gmr, 0.95, fac=1.10)
    t_lim = growth_lim(numpy.ma.log10(t), 0.99, fac=1.05)
    dt_lim = growth_lim(numpy.ma.log10(dt), 0.99, fac=1.05)
    tn_lim = growth_lim(numpy.ma.log10(tn), 0.99, fac=1.05)
    dtn_lim = growth_lim(numpy.ma.log10(dtn), 0.99, fac=1.05)

    sv_lim = growth_lim(svel, 0.90, fac=1.2, midpoint=0)
    so_lim = numpy.power(10., growth_lim(numpy.ma.log10(ssigo), 0.90, fac=1.10))
    sc_lim = growth_lim(ssigcor, 0.90, fac=1.10)
    ss_lim = numpy.power(10., growth_lim(numpy.ma.log10(ssigc), 0.90, fac=1.10))

    so_lim = [ min(so_lim[0],ss_lim[0]), max(so_lim[1],ss_lim[1]) ]
    ss_lim = so_lim.copy()

#    so_lim = [10,400]
#    ss_lim = [10,400]

    # ------------------------------------------------------------------
    # Layout
    left = 0.035
    bott = 0.08
    imwd = 0.18
    hbuf = 0.06
    cbuf = 0.005
    cbwd = 0.01
    vbuf = 0.03
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # S/N
    i,j=0,3
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
#    print(numpy.sum(snr.mask), numpy.prod(snr.shape), numpy.sum(snr.mask)/numpy.prod(snr.shape))
    masked_imshow(fig, ax, cax, snr, extent=extent,
                  norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]), cmap='viridis',
                  zorder=3, cbformat='%.1f', subs=(1.,2.,4.))
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'S/N', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # 68th percentile of fractional residuals
    i,j=1,3
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    print(numpy.sum(r68.mask), numpy.prod(r68.shape), numpy.sum(r68.mask)/numpy.prod(r68.shape))
    masked_imshow(fig, ax, cax, r68, extent=extent,
                  norm=colors.LogNorm(vmin=r68_lim[0], vmax=r68_lim[1]), cmap='viridis_r',
                  zorder=3, cbformat='%.2f', subs=(1.,2.,4.))
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'68% Frac. Resid', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

#    ax.text(0.5, 1.2, '{0}-{1}'.format(plt, ifu), horizontalalignment='center',
#            verticalalignment='center', transform=ax.transAxes, fontsize=12)
#    ax.text(0.5, 1.08, daptype, horizontalalignment='center',
#            verticalalignment='center', transform=ax.transAxes, fontsize=12)

    # ------------------------------------------------------------------
    # 99th percentile of fractional residuals
    i,j=2,3
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    print(numpy.sum(r99.mask), numpy.prod(r99.shape), numpy.sum(r99.mask)/numpy.prod(r99.shape))
    masked_imshow(fig, ax, cax, r99, extent=extent,
                  norm=colors.LogNorm(vmin=r99_lim[0], vmax=r99_lim[1]), cmap='viridis_r',
                  zorder=3, cbformat='%.2f', subs=(1.,2.,4.))
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'99% Frac. Resid', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Reduced chi-square
    i,j=3,3
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    print(numpy.sum(rchi2.mask), numpy.prod(rchi2.shape),
#          numpy.sum(rchi2.mask)/numpy.prod(rchi2.shape))
    masked_imshow(fig, ax, cax, rchi2, extent=extent, vmin=chi_lim[0], vmax=chi_lim[1],
                  cmap='viridis_r', zorder=3, cbformat='%.2f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\chi^2_\nu$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)


    # ------------------------------------------------------------------
    # Signal
    i,j=0,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
#    print(numpy.sum(signal.mask), numpy.prod(signal.shape),
#          numpy.sum(signal.mask)/numpy.prod(signal.shape))
    masked_imshow(fig, ax, cax, signal, extent=extent,
                  norm=colors.LogNorm(vmin=s_lim[0], vmax=s_lim[1]), cmap='viridis',
                  zorder=3, cbformat='%.1f', subs=(1.,2.,4.))
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$I_g$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # A
    i,j=1,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    lga = numpy.ma.log10(a)
#    print(numpy.sum(lga.mask), numpy.prod(lga.shape), numpy.sum(lga.mask)/numpy.prod(lga.shape))
    masked_imshow(fig, ax, cax, lga, extent=extent, vmin=a_lim[0], vmax=a_lim[1], cmap='viridis',
                  zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log A$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Delta A
    i,j=2,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    lgan = numpy.ma.log10(an)
#    print(numpy.sum(lgan.mask), numpy.prod(lgan.shape), numpy.sum(lgan.mask)/numpy.prod(lgan.shape))
    masked_imshow(fig, ax, cax, lgan, extent=extent, vmin=an_lim[0], vmax=an_lim[1],
                  cmap='viridis_r', zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log A_n$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Delta A/Ig
    i,j=3,2
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    lgdan = numpy.ma.log10(dan)
#    print(numpy.sum(lgdan.mask), numpy.prod(lgdan.shape),
#          numpy.sum(lgdan.mask)/numpy.prod(lgdan.shape))
    masked_imshow(fig, ax, cax, lgdan, extent=extent, vmin=dan_lim[0], vmax=dan_lim[1],
                  cmap='viridis_r', zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log(\delta A_n)$', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    # ------------------------------------------------------------------
    # g-r Color
    i,j=0,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
#    print(numpy.sum(gmr.mask), numpy.prod(gmr.shape), numpy.sum(gmr.mask)/numpy.prod(gmr.shape))
    masked_imshow(fig, ax, cax, gmr, extent=extent, vmin=gmr_lim[0], vmax=gmr_lim[1],
                  cmap='RdBu_r', zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$g-r$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # T
    i,j=1,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    lgt = numpy.ma.log10(t)
#    print(numpy.sum(lgt.mask), numpy.prod(lgt.shape), numpy.sum(lgt.mask)/numpy.prod(lgt.shape))
    masked_imshow(fig, ax, cax, lgt, extent=extent, vmin=t_lim[0], vmax=t_lim[1],
                  cmap='viridis', zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log T$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Delta T
    i,j=2,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    lgtn = numpy.ma.log10(tn)
#    print(numpy.sum(lgtn.mask), numpy.prod(lgtn.shape), numpy.sum(lgtn.mask)/numpy.prod(lgtn.shape))
    masked_imshow(fig, ax, cax, lgtn, extent=extent, vmin=tn_lim[0], vmax=tn_lim[1],
                  cmap='viridis_r', zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log T_n$', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Delta T/Ig
    i,j=3,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    lgdtn = numpy.ma.log10(dtn)
#    print(numpy.sum(lgdtn.mask), numpy.prod(lgdtn.shape),
#            numpy.sum(lgdtn.mask)/numpy.prod(lgdtn.shape))
    masked_imshow(fig, ax, cax, lgdtn, extent=extent, vmin=dtn_lim[0], vmax=dtn_lim[1],
                  cmap='viridis_r', zorder=3, cbformat='%.1f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log \delta T_n$', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Stellar velocity
    i,j=0,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
#    print(numpy.sum(svel.mask), numpy.prod(svel.shape), numpy.sum(svel.mask)/numpy.prod(svel.shape))
    masked_imshow(fig, ax, cax, svel, extent=extent, vmin=sv_lim[0], vmax=sv_lim[1],
                  cmap='RdBu_r', zorder=3, cbformat='%.0f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$V_\ast$ (km/s)', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Stellar sigma
    i,j=1,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    print(numpy.sum(ssigo.mask), numpy.prod(ssigo.shape),
#            numpy.sum(ssigo.mask)/numpy.prod(ssigo.shape))
    masked_imshow(fig, ax, cax, ssigo, extent=extent,
                  norm=colors.LogNorm(vmin=so_lim[0], vmax=so_lim[1]),
                  cmap='viridis', zorder=3, cbformat='%.0f', subs=(1.,2.,4.))
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\sigma_{\rm obs}$ (km/s)', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Sigma correction
    i,j=2,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    print(numpy.sum(ssigcor.mask), numpy.prod(ssigcor.shape),
#            numpy.sum(ssigcor.mask)/numpy.prod(ssigcor.shape))
    masked_imshow(fig, ax, cax, ssigcor, extent=extent, vmin=sc_lim[0], vmax=sc_lim[1],
                  cmap='viridis', zorder=3, cbformat='%.0f')
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\sigma_{\rm corr}$ (km/s)', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    # ------------------------------------------------------------------
    # Corrected stellar sigma 
    i,j=3,0
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
    cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
#    print(numpy.sum(ssigc.mask), numpy.prod(ssigc.shape),
#            numpy.sum(ssigc.mask)/numpy.prod(ssigc.shape))
    masked_imshow(fig, ax, cax, ssigo, extent=extent,
                  norm=colors.LogNorm(vmin=ss_lim[0], vmax=ss_lim[1]),
                  cmap='viridis', zorder=3, cbformat='%.0f', subs=(1.,2.,4.))
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.add_patch(map_beam_patch(extent, ax, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\sigma_\ast$ (km/s)', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


def gmr_data(plt, ifu, drpver, redux_path):
    # Get the g-r map from the data cube

    drp_cube_file = os.path.join(*MaNGADataCube.default_paths(plt, ifu, drpver=drpver,
                                                              redux_path=redux_path))
    if not os.path.isfile(drp_cube_file):
        raise FileNotFoundError('{0} does not exist!'.format(drp_cube_file))

    with fits.open(drp_cube_file) as hdu:
        return -2.5*numpy.ma.log10(numpy.ma.MaskedArray(hdu['GIMG'].data,
                                                        mask=numpy.invert(hdu['GIMG'].data>0))
                                    / numpy.ma.MaskedArray(hdu['RIMG'].data,
                                                           mask=numpy.invert(hdu['RIMG'].data>0)))


def rdxqa_data(plt, ifu, plan, drpver, dapver, analysis_path):
    # Get the surface brightness and S/N maps from the
    # ReductionAssessments object
    rdxqa_file = os.path.join(*ReductionAssessment.default_paths(plt, ifu, plan['drpqa_key'],
                                                                 drpver=drpver, dapver=dapver,
                                                                 analysis_path=analysis_path))
    if not os.path.isfile(rdxqa_file):
        raise FileNotFoundError('{0} does not exist!'.format(rdxqa_file))

    with fits.open(rdxqa_file) as hdu:
        spatial_shape = (int(numpy.sqrt(hdu['SPECTRUM'].data['SNR'].size)),)*2
        fgood_map = hdu['SPECTRUM'].data['FGOODPIX'].reshape(spatial_shape).T
        signal_map = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['SIGNAL'].reshape(spatial_shape).T,
                                          mask=numpy.invert(fgood_map > 0))
        snr_map = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['SNR'].reshape(spatial_shape).T,
                                       mask=numpy.invert(fgood_map > 0))
    return signal_map, snr_map

    
def continuum_component_data(plt, ifu, plan, drpver, dapver, analysis_path, signal_map=None,
                             tpl_flux_renorm=None):
    # Get the coefficient data from the StellarContinuumModel object
    sc_file = os.path.join(*StellarContinuumModel.default_paths(plt, ifu, plan['drpqa_key'],
                                                                plan['bin_key'],
                                                                plan['continuum_key'],
                                                                drpver=drpver, dapver=dapver,
                                                                analysis_path=analysis_path))
    hdu = fits.open(sc_file)
    # This only has the bins with stellar continuum measurements
    binid_map = hdu['BINID'].data.copy()

    if signal_map is None:
        signal = numpy.ones(hdu['PAR'].data['BINID_INDEX'].size, dtype=float)
    else:
        # Find the unique indices in the map with valid bins
        unique, unique_indx = numpy.unique(binid_map.ravel(), return_index=True)
        if unique[0] == -1:
            unique_indx = unique_indx[1:]
        # .. and pull these out of the signal map
        signal = signal_map.ravel()[unique_indx]

#    print(hdu['PAR'].data['TPLWGT'].shape)
    w = hdu['PAR'].data['TPLWGT'].copy()
    if tpl_flux_renorm is not None:
        w /= tpl_flux_renorm[None,:]
    wn = w/signal[:,None]
#    t = numpy.sqrt(numpy.sum(numpy.square(w), axis=1))
    t = numpy.sum(numpy.absolute(w), axis=1)
    tn = numpy.sum(numpy.absolute(wn), axis=1)
#    dtw = numpy.median(w / numpy.sum(w, axis=1)[:,None], axis=0)
    dtw = numpy.ones(w.shape[1], dtype=float)
#    dtnw = numpy.median(wn / numpy.sum(wn, axis=1)[:,None], axis=0)
    dtnw = numpy.ones(w.shape[1], dtype=float)
    dt = numpy.ma.sqrt(numpy.sum(dtw[None,:]*numpy.square(w - numpy.ma.median(w, axis=0)[None,:]),
                            axis=1) / numpy.sum(dtw)).filled(0.0)
    dtn = numpy.ma.sqrt(numpy.sum(dtnw[None,:]*numpy.square(wn - numpy.ma.median(wn, axis=0)[None,:]),   
                            axis=1) / numpy.sum(dtnw)).filled(0.0)

    w = hdu['PAR'].data['ADDCOEF'].copy()
    wn = w/signal[:,None]
#    a = numpy.sqrt(numpy.sum(numpy.square(w), axis=1))
    a = numpy.sum(numpy.absolute(w), axis=1)
    an = numpy.sum(numpy.absolute(wn), axis=1)
#    daw = numpy.median(numpy.absolute(w)
#                        / numpy.sum(numpy.absolute(w), axis=1)[:,None], axis=0)
    daw = numpy.ones(w.shape[1], dtype=float)
#    danw = numpy.median(numpy.absolute(wn)
#                        / numpy.sum(numpy.absolute(wn), axis=1)[:,None], axis=0)
    danw = numpy.ones(w.shape[1], dtype=float)
    da = numpy.ma.sqrt(numpy.sum(daw[None,:]*numpy.square(w - numpy.ma.median(w, axis=0)[None,:]),
                            axis=1) / numpy.sum(daw)).filled(0.0)
    dan = numpy.ma.sqrt(numpy.sum(danw[None,:]*numpy.square(wn - numpy.ma.median(wn, axis=0)[None,:]),
                            axis=1) / numpy.sum(danw)).filled(0.0)

    hdu.close()
    del hdu

    return binid_map, t, dt, tn, dtn, a, da, an, dan

    
def maps_data(plt, ifu, plan, drpver, dapver, analysis_path):
    bin_method = SpatiallyBinnedSpectra.define_method(plan['bin_key'])
    sc_method = StellarContinuumModel.define_method(plan['continuum_key'])
    el_method = EmissionLineModel.define_method(plan['elfit_key'])
    method = defaults.dap_method(bin_method['key'], sc_method['fitpar']['template_library_key'],
                                 el_method['continuum_tpl_key'])
    directory_path = defaults.dap_method_path(method, plate=plt, ifudesign=ifu, drpver=drpver,
                                              dapver=dapver, analysis_path=analysis_path)
    maps_file = defaults.dap_file_name(plt, ifu, method, mode='MAPS')

    with fits.open(os.path.join(directory_path, maps_file)) as hdu:

        mask = hdu['BINID'].data[1,:,:] < 0
        r68_map = numpy.ma.MaskedArray(hdu['STELLAR_FOM'].data[3,:,:], mask=mask)
        r99_map = numpy.ma.MaskedArray(hdu['STELLAR_FOM'].data[4,:,:], mask=mask)
        rchi2_map = numpy.ma.MaskedArray(hdu['STELLAR_FOM'].data[2,:,:], mask=mask)
        svel_map = numpy.ma.MaskedArray(hdu['STELLAR_VEL'].data[:,:], mask=mask)
        ssigo_map = numpy.ma.MaskedArray(hdu['STELLAR_SIGMA'].data[:,:], mask=mask)
        ssigcor_map = numpy.ma.MaskedArray(hdu['STELLAR_SIGMACORR'].data[1,:,:], mask=mask)
        ssigc_map = numpy.ma.sqrt( numpy.square(ssigo_map) - numpy.square(ssigcor_map) )

        extent = map_extent(hdu, 'SPX_MFLUX')

    return r68_map, r99_map, rchi2_map, svel_map, ssigo_map, ssigcor_map, ssigc_map, extent


def ppxffit_qa_plot(plt, ifu, plan, drpver=None, redux_path=None, dapver=None, analysis_path=None,
                    tpl_flux_renorm=None):

    # Set the redux and DRP directory paths
    _drpver = defaults.drp_version() if drpver is None else drpver
    _redux_path = defaults.drp_redux_path(drpver=_drpver) if redux_path is None else redux_path
    if not os.path.isdir(_redux_path):
        raise NotADirectoryError('{0} is not a directory.'.format(_redux_path))

    # Set the analysis path
    _dapver = defaults.dap_version() if dapver is None else dapver
    _analysis_path = defaults.dap_analysis_path(drpver=_drpver, dapver=_dapver) \
                            if analysis_path is None else analysis_path
    if not os.path.isdir(_analysis_path):
        raise NotADirectoryError('{0} is not a directory.'.format(_analysis_path))

    # Get the method, it's root directory, and the qa directories
    bin_method = SpatiallyBinnedSpectra.define_method(plan['bin_key'])
    sc_method = StellarContinuumModel.define_method(plan['continuum_key'])
    el_method = EmissionLineModel.define_method(plan['elfit_key'])
    method = defaults.dap_method(bin_method['key'], sc_method['fitpar']['template_library_key'],
                                 el_method['continuum_tpl_key'])
    method_dir = defaults.dap_method_path(method, plate=plt, ifudesign=ifu, drpver=_drpver,
                                          dapver=_dapver, analysis_path=_analysis_path)
    method_qa_dir = defaults.dap_method_path(method, plate=plt, ifudesign=ifu, qa=True,
                                             drpver=_drpver, dapver=_dapver,
                                             analysis_path=_analysis_path)

    # Check that the paths exists
    if not os.path.isdir(method_dir):
        raise NotADirectoryError('{0} is not a directory; run the DAP first!'.format(method_dir))
    if not os.path.isdir(method_qa_dir):
        os.makedirs(method_qa_dir)

    # Get the name of the output file
    ofile = os.path.join(method_qa_dir,
                         'manga-{0}-{1}-MAPS-{2}-ppxffit.png'.format(plt, ifu, method))

    # Get the g-r map
    gmr_map = gmr_data(plt, ifu, _drpver, _redux_path)

    # Get the reduction assessment data
    signal_map, snr_map = rdxqa_data(plt, ifu, plan, _drpver, _dapver, _analysis_path)

    # Get the template weights and additive-polynomial coefficients
    binid_map, t, dt, tn, dtn, a, da, an, dan \
            = continuum_component_data(plt, ifu, plan, _drpver, _dapver, _analysis_path,
                                       signal_map=signal_map, tpl_flux_renorm=tpl_flux_renorm)

#    print('t:', numpy.sum(t > 0), numpy.prod(t.shape))
#    print('dt:', numpy.sum(dt > 0), numpy.prod(dt.shape))
#    print('tn:', numpy.sum(tn > 0), numpy.prod(tn.shape))
#    print('dtn:', numpy.sum(dtn > 0), numpy.prod(dtn.shape))

#    print('a:', numpy.sum(a > 0), numpy.prod(a.shape))
#    print('da:', numpy.sum(da > 0), numpy.prod(da.shape))
#    print('an:', numpy.sum(an > 0), numpy.prod(an.shape))
#    print('dan:', numpy.sum(dan > 0), numpy.prod(dan.shape))

    # Map the continuum component data
    t_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), t, quiet=True)
    t_map = numpy.ma.MaskedArray(t_map, mask=numpy.invert(t_map > 0))

    dt_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), dt, quiet=True)
    dt_map = numpy.ma.MaskedArray(dt_map, mask=numpy.invert(dt_map > 0)) #/t_map

    tn_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), tn, quiet=True)
    tn_map = numpy.ma.MaskedArray(tn_map, mask=numpy.invert(tn_map > 0))

    dtn_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), dtn, quiet=True)
    dtn_map = numpy.ma.MaskedArray(dtn_map, mask=numpy.invert(dtn_map > 0)) #/tn_map

    a_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), a, quiet=True)
    a_map = numpy.ma.MaskedArray(a_map, mask=numpy.invert(a_map > 0))

    da_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), da, quiet=True)
    da_map = numpy.ma.MaskedArray(da_map, mask=numpy.invert(da_map > 0)) #/a_map

    an_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), an, quiet=True)
    an_map = numpy.ma.MaskedArray(an_map, mask=numpy.invert(an_map > 0))

    dan_map = DAPFitsUtil.reconstruct_map(binid_map.shape, binid_map.ravel(), dan, quiet=True)
    dan_map = numpy.ma.MaskedArray(dan_map, mask=numpy.invert(dan_map > 0)) #/an_map

    # Get the remaining output maps from the MAPS file
    r68_map, r99_map, rchi2_map, svel_map, ssigo_map, ssigcor_map, ssigc_map, extent \
            = maps_data(plt, ifu, plan, _drpver, _dapver, _analysis_path)

#    if not os.path.isfile(ofile):
    stellar_continuum_maps(plt, ifu, method, snr_map, r68_map, r99_map, rchi2_map, signal_map,
                           a_map, da_map, an_map, dan_map, gmr_map, t_map, dt_map, tn_map, dtn_map,
                           svel_map, ssigo_map, ssigcor_map, ssigc_map, extent=extent, ofile=ofile)


class PpxfFitQA(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Construct QA plots for pPXF fit', width=width)

        parser.add_argument('plate', type=int, help='plate ID to process')
        parser.add_argument('ifudesign', type=int, help='IFU design to process')

        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('--redux_path', type=str, help='main DRP output path', default=None)
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument('--analysis_path', type=str, help='main DAP output path', default=None)

        parser.add_argument('--plan_file', type=str, help='parameter file with the MaNGA DAP '
                            'execution plan to use instead of the default' , default=None)

        parser.add_argument('--normal_backend', dest='bgagg', action='store_false', default=True)

        parser.add_argument('--template_flux_file', type=str, default=None,
                            help='template renormalization flux file.  Will attempt to read '
                                 'default if not provided.  If no file is provided and the '
                                 'default file does not exist, no renormalization of the '
                                 'templates is performed.')
        return parser

    @staticmethod
    def main(args):
        t = time.perf_counter()

        if args.bgagg:
            pyplot.switch_backend('agg')

        # Get the DAP method types to plot
        analysisplan = AnalysisPlanSet.default() if args.plan_file is None \
                            else AnalysisPlanSet.from_par_file(args.plan_file)

        # Construct the plot for each analysis plan
        for plan in analysisplan:

            # Get the template library keyword
            tpl_renorm_file = args.template_flux_file
            if tpl_renorm_file is None:
                sc_method = StellarContinuumModel.define_method(plan['continuum_key'])
                tpl_key = sc_method['fitpar']['template_library_key']
                library = TemplateLibrary.define_library(tpl_key)
                library_root = os.path.split(library['file_search'])[0]
                tpl_renorm_file = os.path.join(library_root,
                                               '{0}_fluxes.db'.format(tpl_key.lower()))
                if not os.path.isfile(tpl_renorm_file):
                    warnings.warn('Could not find file: {0}'.format(tpl_renorm_file))
                    tpl_renorm_file = None

            tpl_flux_renorm = None if tpl_renorm_file is None \
                                else numpy.genfromtxt(tpl_renorm_file, dtype=float)[:,2]

            # Construct the plot
            ppxffit_qa_plot(args.plate, args.ifudesign, plan, drpver=args.drpver,
                            redux_path=args.redux_path, dapver=args.dapver,
                            analysis_path=args.analysis_path, tpl_flux_renorm=tpl_flux_renorm)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



