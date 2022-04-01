import time

from IPython import embed

import numpy

from matplotlib import pyplot, rc, colors, colorbar, ticker, cm

from astropy.io import fits

from mangadap import dapfits
from mangadap.datacube import DataCube
from mangadap.config.analysisplan import AnalysisPlan
from mangadap.util.fileio import channel_dictionary
from mangadap.proc.util import growth_lim
from mangadap.util.mapping import map_extent, map_beam_patch
from mangadap.util.pkg import load_object

from mangadap.scripts import scriptbase


def init_ax(fig, pos):
    ax = fig.add_axes(pos, facecolor='0.9')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6, direction='in')
    ax.tick_params(which='minor', length=3, direction='in')
    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
#    ax.grid(True, which='minor', color='0.7', zorder=0, linestyle=':')
    return ax


def masked_imshow(fig, ax, cax, data, extent=None, norm=None, vmin=None, vmax=None, cmap='viridis',
                  zorder=0, cbformat=None, subs=None, orientation='vertical', contour_data=None,
                  levels=None, cblocator=None, cbticks=None):
    if numpy.sum(data.mask) != numpy.prod(data.shape):
        if norm is None:
            img = ax.imshow(data, origin='lower', interpolation='nearest', extent=extent,
                            vmin=vmin, vmax=vmax, cmap=cmap, zorder=zorder)
        else:
            img = ax.imshow(data, origin='lower', interpolation='nearest', extent=extent,
                            norm=norm, cmap=cmap, zorder=zorder)

        cb = fig.colorbar(img, cax=cax, orientation=orientation) if cbformat is None else \
                    fig.colorbar(img, cax=cax, format=cbformat, orientation=orientation)
        if subs is not None:
            cb.locator = ticker.LogLocator(base=10, subs=(1.,2.,4.,))
            cb.update_ticks()
        elif cblocator is not None:
            cb.locator = cblocator
            cb.update_ticks()
        elif cbticks is not None:
            cb.set_ticks(cbticks)
            cb.update_ticks()
    else:
        _norm = colors.Normalize(vmin=vmin, vmax=vmax) if norm is None else norm
        cb = colorbar.ColorbarBase(cax, cmap=cm.get_cmap(cmap), norm=_norm)

    if contour_data is not None and levels is not None \
            and numpy.sum(contour_data.mask) != numpy.prod(contour_data.shape):
        cnt = ax.contour(contour_data, origin='lower', extent=extent, colors='0.5',
                         levels=levels, linewidths=1.5, zorder=5)


# Test images:
#  - SPX_MFLUX, SPX_SNR, BINID, BIN_SNR
#  - STELLAR_CONT_FRESID, STELLAR_CONT_CHI2, STELLAR_VEL, STELLAR_SIGMA
#  - EMLINE_SFLUX, EMLINE_GFLUX, EMLINE_GVEL, EMLINE_GSIGMA - all H-alpha
#  - EMLINE_SFLUX (H-beta), EMLINE_GFLUX (H-beta), D4000, Dn4000
def spotcheck_images(cube, plan, method_dir, ref_dir, qa_dir, fwhm=None):
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
        fwhm (:obj:`float`, optional):
            The FWHM in arcsec of the beam.  If None, no beam included in the
            plot.
    """
    # Check that the paths exists
    if not method_dir.exists():
        raise NotADirectoryError(f'{method_dir} does not exist; run the DAP first!')
    if not ref_dir.exists():
        raise NotADirectoryError(f'{ref_dir} does not exist; run the DAP first!')
    if not qa_dir.exists():
        qa_dir.mkdir(parents=True)

    # Get the name of the output file
    ofile = qa_dir / f'{cube.output_root}-MAPS-{plan["key"]}-spotcheck.png'

    bm = dapfits.DAPMapsBitMask()

    _, directory, file = dapfits.default_paths(cube, method=plan['key'], output_path=method_dir)
    ifile = directory / file
    if not ifile.exists():
        raise FileNotFoundError(f'No file: {ifile}')

    print(f'Reading: {ifile}')
    hdu = fits.open(str(ifile))

    # Check everything is finite
    num_ext = len(hdu)
    for i in range(num_ext):
        if hdu[i].data is None:
            continue
        if not numpy.all(numpy.isfinite(hdu[i].data)):
            raise ValueError(f'HDU {hdu[i].name} contains infs or NaNs!')

    # Build the column dictionaries
    emline = channel_dictionary(hdu, 'EMLINE_GFLUX')
    specindex = channel_dictionary(hdu, 'SPECINDEX')

    extent = map_extent(hdu, 'SPX_MFLUX')

    # Get the data to plot
    spxflx = numpy.ma.MaskedArray(hdu['SPX_MFLUX'].data.copy())
    spxsnr = numpy.ma.MaskedArray(hdu['SPX_SNR'].data.copy())
    binid = numpy.ma.MaskedArray(hdu['BINID'].data.copy(), mask=hdu['BINID'].data<0)
    binsnr = numpy.ma.MaskedArray(hdu['BIN_SNR'].data.copy(), mask=hdu['BIN_MFLUX_MASK'].data>0)

    # 68% growth of the absolute value of the fractional residuals
    scfres = numpy.ma.MaskedArray(hdu['STELLAR_FOM'].data.copy()[3,:,:],
                                  mask=numpy.invert(hdu['STELLAR_FOM'].data[3,:,:] > 0))
    # Reduced chi-square
    scrchi = numpy.ma.MaskedArray(hdu['STELLAR_FOM'].data.copy()[2,:,:],
                                  mask=numpy.invert(hdu['STELLAR_FOM'].data[2,:,:] > 0))

    strvel = numpy.ma.MaskedArray(hdu['STELLAR_VEL'].data.copy(),
                                  mask=bm.flagged(hdu['STELLAR_VEL_MASK'].data, 'DONOTUSE'))
    ustrvel = numpy.ma.MaskedArray(hdu['STELLAR_VEL'].data.copy(),
                                  mask=numpy.invert(bm.flagged(hdu['STELLAR_VEL_MASK'].data,
                                                               'UNRELIABLE')))
    ustrvel[numpy.invert(numpy.ma.getmaskarray(ustrvel))] = 0.0
    strsig = numpy.ma.MaskedArray(hdu['STELLAR_SIGMA'].data.copy(),
                                  mask=bm.flagged(hdu['STELLAR_SIGMA_MASK'].data, 'DONOTUSE'))
    ustrsig = numpy.ma.MaskedArray(hdu['STELLAR_SIGMA'].data.copy(),
                                  mask=numpy.invert(bm.flagged(hdu['STELLAR_SIGMA_MASK'].data,
                                                               'UNRELIABLE')))
    ustrsig[numpy.invert(numpy.ma.getmaskarray(ustrsig))] = 0.0

    hasflx = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX'].data.copy()[emline['Ha-6564'],:,:],
            mask=bm.flagged(hdu['EMLINE_SFLUX_MASK'].data[emline['Ha-6564'],:,:], 'DONOTUSE'))
    uhasflx = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX'].data.copy()[emline['Ha-6564'],:,:],
          mask=numpy.invert(bm.flagged(hdu['EMLINE_SFLUX_MASK'].data[emline['Ha-6564'],:,:],
                                       'UNRELIABLE')))
    uhasflx[numpy.invert(numpy.ma.getmaskarray(uhasflx))] = 0.0
    hagflx = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data.copy()[emline['Ha-6564'],:,:],
          mask=bm.flagged(hdu['EMLINE_GFLUX_MASK'].data[emline['Ha-6564'],:,:], 'DONOTUSE'))
    uhagflx = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data.copy()[emline['Ha-6564'],:,:],
          mask=numpy.invert(bm.flagged(hdu['EMLINE_GFLUX_MASK'].data[emline['Ha-6564'],:,:],
                                       'UNRELIABLE')))
    uhagflx[numpy.invert(numpy.ma.getmaskarray(uhagflx))] = 0.0
    hagvel = numpy.ma.MaskedArray(hdu['EMLINE_GVEL'].data.copy()[emline['Ha-6564'],:,:],
          mask=bm.flagged(hdu['EMLINE_GVEL_MASK'].data[emline['Ha-6564'],:,:], 'DONOTUSE'))
    uhagvel = numpy.ma.MaskedArray(hdu['EMLINE_GVEL'].data.copy()[emline['Ha-6564'],:,:],
          mask=numpy.invert(bm.flagged(hdu['EMLINE_GVEL_MASK'].data[emline['Ha-6564'],:,:],
                                       'UNRELIABLE')))
    uhagvel[numpy.invert(numpy.ma.getmaskarray(uhagvel))] = 0.0
    hagsig = numpy.ma.MaskedArray(hdu['EMLINE_GSIGMA'].data.copy()[emline['Ha-6564'],:,:],
        mask=bm.flagged(hdu['EMLINE_GSIGMA_MASK'].data[emline['Ha-6564'],:,:], 'DONOTUSE'))
    uhagsig = numpy.ma.MaskedArray(hdu['EMLINE_GSIGMA'].data.copy()[emline['Ha-6564'],:,:],
        mask=numpy.invert(bm.flagged(hdu['EMLINE_GSIGMA_MASK'].data[emline['Ha-6564'],:,:],
                                     'UNRELIABLE')))
    uhagsig[numpy.invert(numpy.ma.getmaskarray(uhagsig))] = 0.0

    hbsflx = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX'].data.copy()[emline['Hb-4862'],:,:],
          mask=bm.flagged(hdu['EMLINE_SFLUX_MASK'].data[emline['Hb-4862'],:,:], 'DONOTUSE'))
    uhbsflx = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX'].data.copy()[emline['Hb-4862'],:,:],
          mask=numpy.invert(bm.flagged(hdu['EMLINE_SFLUX_MASK'].data[emline['Hb-4862'],:,:],
                                       'UNRELIABLE')))
    uhbsflx[numpy.invert(numpy.ma.getmaskarray(uhbsflx))] = 0.0
    hbgflx = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data.copy()[emline['Hb-4862'],:,:],
          mask=bm.flagged(hdu['EMLINE_GFLUX_MASK'].data[emline['Hb-4862'],:,:], 'DONOTUSE'))
    uhbgflx = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data.copy()[emline['Hb-4862'],:,:],
          mask=numpy.invert(bm.flagged(hdu['EMLINE_GFLUX_MASK'].data[emline['Hb-4862'],:,:],
                                       'UNRELIABLE')))
    uhbgflx[numpy.invert(numpy.ma.getmaskarray(uhbgflx))] = 0.0

    specindex_available = True
    try:
        if hdu['SPECINDEX'].data is None:
            pass
    except:
        warnings.warn('No spectral indices in maps file.')
        specindex_available = False
    if specindex_available:
        d4000 = numpy.ma.MaskedArray(hdu['SPECINDEX'].data.copy()[specindex['D4000'],:,:],
                                mask=bm.flagged(hdu['SPECINDEX_MASK'].data[specindex['D4000'],:,:],
                                                'DONOTUSE'))
        ud4000 = numpy.ma.MaskedArray(hdu['SPECINDEX'].data.copy()[specindex['D4000'],:,:],
                                mask=numpy.invert(bm.flagged(hdu['SPECINDEX_MASK'].data[
                                                      specindex['D4000'],:,:], 'UNRELIABLE')))
        ud4000[numpy.invert(numpy.ma.getmaskarray(ud4000))] = 0.0
        dn4000 = numpy.ma.MaskedArray(hdu['SPECINDEX'].data.copy()[specindex['Dn4000'],:,:],
                                mask=bm.flagged(hdu['SPECINDEX_MASK'].data[specindex['Dn4000'],:,:],
                                                'DONOTUSE'))
        udn4000 = numpy.ma.MaskedArray(hdu['SPECINDEX'].data.copy()[specindex['Dn4000'],:,:],
                                mask=numpy.invert(bm.flagged(
                                        hdu['SPECINDEX_MASK'].data[specindex['Dn4000'],:,:],
                                        'UNRELIABLE')))
        udn4000[numpy.invert(numpy.ma.getmaskarray(udn4000))] = 0.0
    else:
        output_shape = spxflx.shape
        d4000 = numpy.ma.masked_all(output_shape)
        ud4000 = numpy.ma.masked_all(output_shape)
        dn4000 = numpy.ma.masked_all(output_shape)
        udn4000 = numpy.ma.masked_all(output_shape)

    # Get the limits to apply
    x_order = 1 if extent[0] < extent[1] else -1
    Dx = x_order*(extent[1] - extent[0])
    y_order = 1 if extent[2] < extent[3] else -1
    Dy = y_order*(extent[3] - extent[2])
    # Force the image panels to be square; assumes extent has the same units in
    # both dimensions.
    D = max(Dx, Dy)

    x_lim = (extent[0]+extent[1])/2 + D * x_order * numpy.array([-1,1])/2
    y_lim = (extent[2]+extent[3])/2 + D * y_order * numpy.array([-1,1])/2

    flux_lim = numpy.power(10, growth_lim(numpy.ma.log10(spxflx), 0.90, fac=1.05))
    t = numpy.ma.append(numpy.ma.log10(spxsnr), numpy.ma.log10(binsnr))
    snr_lim = numpy.power(10, growth_lim(t, 0.90, fac=1.05))
    bin_lim = [ numpy.ma.amin(binid), numpy.ma.amax(binid) ]

    res_lim = numpy.power(10, growth_lim(numpy.ma.log10(scfres), 0.90, fac=1.05))
    chi_lim = growth_lim(scrchi, 0.90, fac=1.05)
    t = numpy.ma.append(strvel, hagvel)
    vel_lim = growth_lim(t, 0.85, fac=1.10, midpoint=0.0)
    if vel_lim[0] < -350:
        vel_lim[0] = -350
    if vel_lim[1] > 350:
        vel_lim[1] = 350
    t = numpy.ma.append(numpy.ma.log10(strsig), numpy.ma.log10(hagsig))
    sig_lim = numpy.power(10, growth_lim(t, 0.90, fac=1.05))

    t = numpy.ma.append( numpy.ma.append(numpy.ma.log10(hasflx), numpy.ma.log10(hagflx)),
                         numpy.ma.append(numpy.ma.log10(hbsflx), numpy.ma.log10(hbgflx)) )
#    t = numpy.ma.log10(t)
    hflx_lim = numpy.power(10, growth_lim(t, 0.95, fac=1.10))
#    hflx_lim = [ -0.3, 2.2 ]
    t = numpy.ma.append(d4000, dn4000)
#    d4000_lim = growth_lim(t, 0.90, 1.05)
    d4000_lim = [ 1.0, 2.5 ]

    font = { 'size' : 6 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    dx = 0.22
    left = 0.05
    top=0.94
    dw = 0.005

    snr_levels = [0.5, 1.0, 1.5, 2.0 ]

    ax = init_ax(fig, [left, top-dx, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx*0.65, top-0.07*dx, 0.27*dx, 0.01])
    cbticks = [flux_lim[0]*1.1, flux_lim[1]/1.1]
    masked_imshow(fig, ax, cax, spxflx, extent=extent,
                  norm=colors.LogNorm(vmin=flux_lim[0], vmax=flux_lim[1]), cmap='YlGnBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter())
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'g-band S (spx)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)

    ax = init_ax(fig, [left+dx+dw, top-dx, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx+dw+dx*0.65, top-0.07*dx, 0.27*dx, 0.01])
    cbticks = [snr_lim[0]*1.1, snr_lim[1]/1.1]
    masked_imshow(fig, ax, cax, spxsnr, extent=extent,
                  norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]), cmap='YlGnBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'g-band S/N (spx)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
   
#    ax.text(1.0, 1.1, r'{0}-{1}; {2}'.format(plate,ifudesign, hdu['PRIMARY'].header['MANGAID']),
#            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
#            fontsize=20)
    ax.text(1.0, 1.1, f'{cube.name}; {plan["key"]}', ha='center', va='center',
            transform=ax.transAxes, fontsize=20)

    ax = init_ax(fig, [left+2*dx+2*dw, top-dx, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+2*dx+2*dw+dx*0.65, top-0.07*dx, 0.27*dx, 0.01])
    masked_imshow(fig, ax, cax, binid[0,:,:], extent=extent,
                  vmin=bin_lim[0], vmax=bin_lim[1], cmap='CMRmap',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbformat=ticker.ScalarFormatter(),
                  cbticks=bin_lim) 
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'bin ID', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
 
    ax = init_ax(fig, [left+3*dx+3*dw, top-dx, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+3*dx+3*dw+dx*0.65, top-0.07*dx, 0.27*dx, 0.01])
    cbticks = [snr_lim[0]*1.1, snr_lim[1]/1.1]
    masked_imshow(fig, ax, cax, binsnr, extent=extent,
                  norm=colors.LogNorm(vmin=snr_lim[0], vmax=snr_lim[1]), cmap='YlGnBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'g-band S/N (bin)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)


    ax = init_ax(fig, [left, top-2*dx-dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx*0.65, top-dx-dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [res_lim[0]*1.1, res_lim[1]/1.1]
    masked_imshow(fig, ax, cax, scfres, extent=extent,
                  norm=colors.LogNorm(vmin=res_lim[0], vmax=res_lim[1]), cmap='inferno',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'pPXF Frac. Resid.', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)

    ax = init_ax(fig, [left+dx+dw,  top-2*dx-dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx+dw+dx*0.65, top-dx-dw-0.07*dx, 0.27*dx, 0.01])
    cbticks=[numpy.round(chi_lim[0]+0.01, decimals=2), numpy.round(chi_lim[1]-0.01, decimals=2)]
    masked_imshow(fig, ax, cax, scrchi, extent=extent,
                  vmin=chi_lim[0], vmax=chi_lim[1], cmap='inferno',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'${\rm pPXF}\ \chi^2_\nu$', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
    
    ax = init_ax(fig, [left+2*dx+2*dw,  top-2*dx-dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+2*dx+2*dw+dx*0.65, top-dx-dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [vel_lim[0]/1.1, 0.0, vel_lim[1]/1.1]
    masked_imshow(fig, ax, cax, strvel, extent=extent,
                  vmin=vel_lim[0], vmax=vel_lim[1], cmap='RdBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'$V_\ast$', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
    
    ax = init_ax(fig, [left+3*dx+3*dw,  top-2*dx-dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+3*dx+3*dw+dx*0.65, top-dx-dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [sig_lim[0]*1.1, numpy.power(10, numpy.mean(numpy.log10(sig_lim))), sig_lim[1]/1.1]
    masked_imshow(fig, ax, cax, strsig, extent=extent,
                  norm=colors.LogNorm(vmin=sig_lim[0], vmax=sig_lim[1]), cmap='viridis',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'$\sigma_\ast$', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)



    ax = init_ax(fig, [left, top-3*dx-2*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx*0.65, top-2*dx-2*dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [hflx_lim[0]*1.1, hflx_lim[1]/1.1]
    masked_imshow(fig, ax, cax, hasflx, extent=extent,
                  norm=colors.LogNorm(vmin=hflx_lim[0], vmax=hflx_lim[1]), cmap='inferno',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'${\rm H}\alpha\ {\rm flux}$ (sum)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)

    ax = init_ax(fig, [left+dx+dw,  top-3*dx-2*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx+dw+dx*0.65, top-2*dx-2*dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [hflx_lim[0]*1.1, hflx_lim[1]/1.1]
    masked_imshow(fig, ax, cax, hagflx, extent=extent,
                  norm=colors.LogNorm(vmin=hflx_lim[0], vmax=hflx_lim[1]), cmap='inferno',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'${\rm H}\alpha\ {\rm flux}$ (Gauss)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
    
    ax = init_ax(fig, [left+2*dx+2*dw,  top-3*dx-2*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+2*dx+2*dw+dx*0.65, top-2*dx-2*dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [vel_lim[0]/1.1, 0.0, vel_lim[1]/1.1]
    masked_imshow(fig, ax, cax, hagvel, extent=extent,
                  vmin=vel_lim[0], vmax=vel_lim[1], cmap='RdBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'$V_{{\rm H}\alpha}$', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
    
    ax = init_ax(fig, [left+3*dx+3*dw,  top-3*dx-2*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+3*dx+3*dw+dx*0.65, top-2*dx-2*dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [sig_lim[0]*1.1, numpy.power(10, numpy.mean(numpy.log10(sig_lim))), sig_lim[1]/1.1]
    masked_imshow(fig, ax, cax, hagsig, extent=extent,
                  norm=colors.LogNorm(vmin=sig_lim[0], vmax=sig_lim[1]), cmap='viridis',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'$\sigma_{{\rm H}\alpha}$', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)



    ax = init_ax(fig, [left, top-4*dx-3*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    cax = fig.add_axes([left+dx*0.65, top-3*dx-3*dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [hflx_lim[0]*1.1, hflx_lim[1]/1.1]
    masked_imshow(fig, ax, cax, hbsflx, extent=extent,
                  norm=colors.LogNorm(vmin=hflx_lim[0], vmax=hflx_lim[1]), cmap='inferno',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'${\rm H}\beta\ {\rm flux}$ (sum)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)

    ax = init_ax(fig, [left+dx+dw,  top-4*dx-3*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+dx+dw+dx*0.65, top-3*dx-3*dw-0.07*dx, 0.27*dx, 0.01])
    cbticks = [hflx_lim[0]*1.1, hflx_lim[1]/1.1]
    masked_imshow(fig, ax, cax, hbgflx, extent=extent,
                  norm=colors.LogNorm(vmin=hflx_lim[0], vmax=hflx_lim[1]), cmap='inferno',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=cbticks, cbformat=ticker.ScalarFormatter()) #'%.0f')
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'${\rm H}\beta\ {\rm flux}$ (Gauss)', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)

    ax = init_ax(fig, [left+2*dx+2*dw,  top-4*dx-3*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+2*dx+2*dw+dx*0.65, top-3*dx-3*dw-0.07*dx, 0.27*dx, 0.01])
    masked_imshow(fig, ax, cax, d4000, extent=extent,
                  vmin=d4000_lim[0], vmax=d4000_lim[1], cmap='RdBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=d4000_lim, cbformat=ticker.ScalarFormatter())
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'D4000', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)
    
    ax = init_ax(fig, [left+3*dx+3*dw,  top-4*dx-3*dw, dx, dx])
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([left+3*dx+3*dw+dx*0.65, top-3*dx-3*dw-0.07*dx, 0.27*dx, 0.01])
    masked_imshow(fig, ax, cax, dn4000, extent=extent,
                  vmin=d4000_lim[0], vmax=d4000_lim[1], cmap='RdBu_r',
                  zorder=4, orientation='horizontal', contour_data=numpy.ma.log10(spxsnr),
                  levels=snr_levels, cbticks=d4000_lim, cbformat=ticker.ScalarFormatter())
    if fwhm is not None:
        ax.add_patch(map_beam_patch(extent, ax, fwhm=fwhm, facecolor='0.7', edgecolor='k', zorder=4))
    ax.text(0.99, 0.05, r'Dn4000', horizontalalignment='right',
            verticalalignment='center', transform=ax.transAxes)

    if ofile is not None:
        fig.canvas.print_figure(ofile, bbox_inches='tight')
        print('Writing: {0}'.format(ofile))
    else:
        pyplot.show()

    fig.clear()
    pyplot.close(fig)
    

class SpotcheckDapMaps(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'spotcheck_dap_maps'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Construct a QA plot to spotcheck DAP results',
                                    width=width)

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

        parser.add_argument('-b', '--beam', type=float, default=None, help='Beam FWHM for plot.')

        parser.add_argument('--normal_backend', dest='bgagg', action='store_false', default=True)

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
            spotcheck_images(cube, plan[key], plan.method_path(plan_index=i),
                             plan.method_path(plan_index=i, ref=True),
                             plan.method_path(plan_index=i, qa=True), fwhm=args.beam)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



