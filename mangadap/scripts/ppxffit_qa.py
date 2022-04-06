import time
import warnings

from IPython import embed

import numpy
from scipy import interpolate
from matplotlib import pyplot, colors, rc, colorbar, ticker, cm

from astropy.io import fits

from mangadap.datacube import DataCube
from mangadap.config.analysisplan import AnalysisPlan
from mangadap.proc.templatelibrary import TemplateLibraryDef
from mangadap.proc.reductionassessments import ReductionAssessment, ReductionAssessmentDef
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectraDef
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelDef
from mangadap.proc.util import growth_lim
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.mapping import map_extent, map_beam_patch
from mangadap.util.pkg import load_object
from mangadap.config import defaults
from mangadap import dapfits

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

# TODO:
#   - Add a buffer keyword that increases the size of the image panels
#   - Use astropy image plotting (with wcs) tools instead?
def stellar_continuum_maps(name, snr, r68, r99, rchi2, signal, a, da, an, dan,
                           gmr, t, dt, tn, dtn, svel, ssigo, ssigcor, ssigc, extent=None,
                           fwhm=2.5, ofile=None):

    font = { 'size' : 6 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(2*w,2*h))

    x_order = 1 if extent[0] < extent[1] else -1
    Dx = x_order*(extent[1] - extent[0])
    y_order = 1 if extent[2] < extent[3] else -1
    Dy = y_order*(extent[3] - extent[2])
    # Force the image panels to be square; assumes extent has the same units in
    # both dimensions.
    D = max(Dx, Dy)

    x_lim = (extent[0]+extent[1])/2 + D * x_order * numpy.array([-1,1])/2
    y_lim = (extent[2]+extent[3])/2 + D * y_order * numpy.array([-1,1])/2

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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'68% Frac. Resid', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

#    ax.text(0.5, 1.2, name, horizontalalignment='center',
#            verticalalignment='center', transform=ax.transAxes, fontsize=12)
#    ax.text(0.5, 1.08, method, horizontalalignment='center',
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\log(\delta A_n)$', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    # ------------------------------------------------------------------
    # g-r Color
    i,j=0,1
    ax = init_image_ax(fig, [left+i*(imwd+hbuf), bott+j*(imwd+vbuf), imwd, imwd ])
#    print(numpy.sum(gmr.mask), numpy.prod(gmr.shape), numpy.sum(gmr.mask)/numpy.prod(gmr.shape))
    if gmr is not None:
        cax = fig.add_axes([left+i*(imwd+hbuf)+imwd+cbuf, bott+j*(imwd+vbuf), cbwd, imwd ])
        masked_imshow(fig, ax, cax, gmr, extent=extent, vmin=gmr_lim[0], vmax=gmr_lim[1],
                      cmap='RdBu_r', zorder=3, cbformat='%.1f')
    else:
        ax.text(0.5, 0.5, 'No g-r data', ha='center', va='center', transform=ax.transAxes)
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
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
    ax.set_xlim(x_lim) #extent[:2])
    ax.set_ylim(y_lim) #extent[2:])
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.1, 0.95, r'$\sigma_\ast$ (km/s)', horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes)

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


# TODO: Should make this a datacube member function
def gmr_data(cube, min_frac=0.8):
    """
    Read or construct a g-r for the provided datacube.
    """
    # First see if there are GIMG and RIMG extensions in the input datacube
    # (requires a MaNGA DRP-like datamodel)
    with fits.open(str(cube.file_path)) as hdu:
        ext_list = [h.name for h in hdu]
        if 'GIMG' in ext_list and 'RIMG' in ext_list:
            return -2.5*numpy.ma.log10(numpy.ma.MaskedArray(hdu['GIMG'].data,
                                                            mask=numpy.invert(hdu['GIMG'].data>0))
                                    / numpy.ma.MaskedArray(hdu['RIMG'].data,
                                                           mask=numpy.invert(hdu['RIMG'].data>0)))

    # Next, try to create it from scratch
    g_file = defaults.dap_data_root() / 'filter_response' / 'gunn_2001_g_response.db'
    if not g_file.exists():
        warnings.warn(f'Could not find file {g_file}; cannot produce rough g-r image.')
        return None

    r_file = defaults.dap_data_root() / 'filter_response' / 'gunn_2001_r_response.db'
    if not r_file.exists():
        warnings.warn(f'Could not find file {r_file}; cannot produce rough g-r image.')
        return None

    flux = cube.copy_to_masked_array()
    # Fix the ordering
    flux = flux.reshape(cube.flux.shape).T.reshape(npix,-1)
    gpm = numpy.logical_not(numpy.ma.getmaskarray(flux)).astype(float)

    # Wavelengths in these files are in angstroms
    g_wave, g_eff = numpy.genfromtxt(str(g_file)).T
    indx = g_eff > 0.
    g_lo = numpy.amin(g_wave[indx])
    g_hi = numpy.amax(g_wave[indx])
    g_wgt = interpolate.interp1d(g_wave[indx], g_eff[indx], bounds_error=False,
                                 fill_value=0.)(cube.wave)[None,:]*gpm
    g_mask = numpy.sum((g_wgt[:, 1:]>0).astype(float) * numpy.diff(cube.wave)[None, :],
                       axis=1).reshape(cube.spatial_shape[::-1]) < min_frac*(g_hi-g_lo)
    if numpy.all(g_mask):
        warnings.warn('Datacube has insufficient coverage of the g-band to construct a g-r image.')
        return None

    r_wave, r_eff = numpy.genfromtxt(str(r_file)).T
    indx = r_eff > 0.
    r_lo = numpy.amin(r_wave[indx])
    r_hi = numpy.amax(r_wave[indx])
    r_wgt = interpolate.interp1d(r_wave[indx], r_eff[indx], bounds_error=False,
                                 fill_value=0.)(cube.wave)[None,:]*gpm

    r_mask = numpy.sum((r_wgt[:, 1:]>0).astype(float) * numpy.diff(cube.wave)[None, :],
                       axis=1).reshape(cube.spatial_shape) < min_frac*(r_hi-r_lo)
    if numpy.all(r_mask):
        warnings.warn('Datacube has insufficient coverage of the r-band to construct a g-r image.')
        return None

    g = numpy.sum(g_wgt[:, 1:] * flux[:, 1:] * numpy.diff(cube.wave)[None,:],
                       axis=1).reshape(cube.spatial_shape[::-1])
    r = numpy.sum(r_wgt[:, 1:] * flux[:, 1:] * numpy.diff(cube.wave)[None,:],
                       axis=1).reshape(cube.spatial_shape[::-1])
    return -2.5*numpy.ma.log10(numpy.ma.divide(numpy.ma.MaskedArray(g, mask=g_mask),
                                               numpy.ma.MaskedArray(r, mask=r_mask)))


def rdxqa_data(cube, rdxqa_method, output_path):
    # Get the surface brightness and S/N maps from the
    # ReductionAssessments object

    directory, file = ReductionAssessment.default_paths(cube, rdxqa_method['key'],
                                                        output_path=output_path)
    rdxqa_file = directory / file
    if not rdxqa_file.exists():
        raise FileNotFoundError(f'{rdxqa_file} does not exist!')

    with fits.open(str(rdxqa_file)) as hdu:
        spatial_shape = (int(numpy.sqrt(hdu['SPECTRUM'].data['SNR'].size)),)*2
        fgood_map = hdu['SPECTRUM'].data['FGOODPIX'].reshape(spatial_shape).T
        signal_map = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['SIGNAL'].reshape(spatial_shape).T,
                                          mask=numpy.invert(fgood_map > 0))
        snr_map = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['SNR'].reshape(spatial_shape).T,
                                       mask=numpy.invert(fgood_map > 0))
    return signal_map, snr_map

    
def continuum_component_data(cube, rdxqa_method, binning_method, continuum_method, output_path,
                             signal_map=None, tpl_flux_renorm=None):

    # Get the coefficient data from the StellarContinuumModel object
    directory, file = StellarContinuumModel.default_paths(cube, continuum_method['key'],
                                                          rdxqa_method['key'],
                                                          binning_method['key'],
                                                          output_path=output_path)
    sc_file = directory / file
    if not sc_file.exists():
        raise FileNotFoundError(f'{sc_file} does not exist!')

    hdu = fits.open(str(sc_file))
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

    
def maps_data(maps_file):
    if not maps_file.exists():
        raise FileNotFoundError(f'{maps_file} does not exist!')

    with fits.open(str(maps_file)) as hdu:

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


def ppxffit_qa_plot(cube, plan, method_dir, ref_dir, qa_dir, fwhm=None, tpl_flux_renorm=None):
    """
    Construct a QA plot for the PPXFFit results.

    Args:
        cube (:class:`~mangadap.datacube.datacube.DataCube`):
            Analyzed data cube
        plan (:obj:`dict`, :class:`~mangadap.par.parset.ParSet`):
            Object with analysis plan keywords.
        method_dir (`Path`_):
            Top-level "method" directory with all the DAP files.
        ref_dir (`Path`_):
            Directory with all the DAP reference files.
        qa_dir (`Path`_):
            Directory for the QA figures
        tpl_flux_renorm (`numpy.ndarray`_, optional):
            Renormalization coefficients for all the templates.  Shape must
            match the number of templates for the stellar conntinuum fit.
    """
    # Check that the paths exists
    if not method_dir.exists():
        raise NotADirectoryError(f'{method_dir} does not exist; run the DAP first!')
    if not ref_dir.exists():
        raise NotADirectoryError(f'{ref_dir} does not exist; run the DAP first!')
    if not qa_dir.exists():
        qa_dir.mkdir(parents=True)

    # Get the name of the output file
    ofile = qa_dir / f'{cube.output_root}-MAPS-{plan["key"]}-ppxffit.png'

    # Get the g-r map
    gmr_map = gmr_data(cube)

    # Get the reduction assessment data
    rdxqa_method = ReductionAssessmentDef.from_dict(plan['rdxqa'])
    signal_map, snr_map = rdxqa_data(cube, rdxqa_method, ref_dir)

    # Get the template weights and additive-polynomial coefficients
    binning_method = SpatiallyBinnedSpectraDef.from_dict(plan['binning'])
    continuum_method = StellarContinuumModelDef.from_dict(plan['continuum'])
    binid_map, t, dt, tn, dtn, a, da, an, dan \
            = continuum_component_data(cube, rdxqa_method, binning_method, continuum_method,
                                       ref_dir, signal_map=signal_map,
                                       tpl_flux_renorm=tpl_flux_renorm)

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

    _, directory, file = dapfits.default_paths(cube, method=plan['key'], output_path=method_dir)
    maps_file = directory / file
    
    # Get the remaining output maps from the MAPS file
    r68_map, r99_map, rchi2_map, svel_map, ssigo_map, ssigcor_map, ssigc_map, extent \
            = maps_data(maps_file)

    stellar_continuum_maps(cube.name, snr_map, r68_map, r99_map, rchi2_map, signal_map,
                           a_map, da_map, an_map, dan_map, gmr_map, t_map, dt_map, tn_map, dtn_map,
                           svel_map, ssigo_map, ssigcor_map, ssigc_map, fwhm=fwhm, extent=extent,
                           ofile=ofile)


class PpxfFitQA(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Construct QA plots for pPXF fit.', width=width)

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-c', '--config', type=str, default=None,
                           help='Configuration file used to instantiate the relevant DataCube '
                                'derived class.')
        group.add_argument('-f', '--cubefile', type=str, default=None,
                           help='Name of the file with the datacube data.  Must be possible to '
                                'instantiate the relevant DataCube derived class directly from '
                                'the file only.')
        parser.add_argument('--cube_module', nargs='*',
                            default='mangadap.datacube.MaNGADataCube',
                            help='The name of the module that contains the DataCube derived '
                                 'class used to read the data.')

        parser.add_argument('-p', '--plan', type=str,
                            help='SDSS parameter file with analysis plan.  If not provided, a '
                                 'default plan is used.')
        parser.add_argument('--plan_module', nargs='*',
                            default='mangadap.config.manga.MaNGAAnalysisPlan',
                            help='The name of the module used to define the analysis plan and '
                                 'the output paths.')

        parser.add_argument('-o', '--output_path', type=str, default=None,
                            help='Top-level directory for the DAP output files; default path is '
                                 'set by the provided analysis plan object (see plan_module).')

        # TODO: For MaNGA, pull fwhm from the header.  Put FWHM in cube.meta?
        parser.add_argument('-b', '--beam', type=float, default=None, help='Beam FWHM for plot.')

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

        # Instantiate the DataCube
        #   - Import the module used to read the datacube
        if isinstance(args.cube_module, list) and len(args.cube_module) > 2:
            raise ValueError('Provided cube module must be one or two strings.')
        if isinstance(args.cube_module, str) or len(args.cube_module) == 1:
            UserDataCube = load_object(args.cube_module if isinstance(args.cube_module, str) 
                                       else args.cube_module[0])
        else:
            UserDataCube = load_object(args.cube_module[0], obj=args.cube_module[1])
        #   - Check that the class is derived from DataCube
        if not issubclass(UserDataCube, DataCube):
            raise TypeError('Defined cube object must subclass from mangadap.datacube.DataCube.')

        #   - Import the module used to set the analysis plan
        if isinstance(args.plan_module, list) and len(args.plan_module) > 2:
            raise ValueError('Provided plan module must be one or two strings.')
        if isinstance(args.plan_module, str) or len(args.plan_module) == 1:
            UserPlan = load_object(args.plan_module if isinstance(args.plan_module, str) 
                                       else args.plan_module[0])
        else:
            UserPlan = load_object(args.plan_module[0], obj=args.plan_module[1])
        #   - Check that the class is derived from AnalysisPlan
        if not issubclass(UserPlan, AnalysisPlan):
            raise TypeError('Defined plan object must subclass from '
                            'mangadap.config.analysisplan.AnalysisPlan')

        #   - Instantiate using either the datacube file directly or a
        #     configuration file
        cube = UserDataCube(args.cubefile) if args.config is None \
                    else UserDataCube.from_config(args.config)

        # Read the analysis plan
        plan = UserPlan.default(cube=cube, analysis_path=args.output_path) if args.plan is None \
                    else UserPlan.from_toml(args.plan, cube=cube, analysis_path=args.output_path)

        # Construct the plot for each analysis plan
        for i, key in enumerate(plan.keys()):

            # Get the template library keyword
            if args.template_flux_file is None:
                sc_method = StellarContinuumModelDef.from_dict(plan[key]['continuum'])
                tpl_key = sc_method['fitpar']['template_library']['key']
                library_root = sc_method["fitpar"]["template_library"].get_file_list()[0].parent
                tpl_renorm_file = library_root / f'{tpl_key.lower()}_fluxes.db'
                #sc_method = StellarContinuumModel.define_method(plan[i]['continuum_key'])
                #tpl_key = sc_method['fitpar']['template_library_key']
                #library = TemplateLibrary.define_library(tpl_key)
                #library_root = library.get_file_list()[0].parent
                #tpl_renorm_file = library_root / f'{tpl_key.lower()}_fluxes.db'
            else:
                tpl_renorm_file = Path(args.template_flux_file).resolve()

            if not tpl_renorm_file.exists():
                warnings.warn('Could not find file with template renormalizations: '
                                f'{tpl_renorm_file}')
                tpl_renorm_file = None

            tpl_flux_renorm = None if tpl_renorm_file is None \
                                else numpy.genfromtxt(tpl_renorm_file, dtype=float)[:,2]

            # Construct the plot
            ppxffit_qa_plot(cube, plan[key], plan.method_path(plan_index=i),
                            plan.method_path(plan_index=i, ref=True),
                            plan.method_path(plan_index=i, qa=True), fwhm=args.beam,
                            tpl_flux_renorm=tpl_flux_renorm)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



