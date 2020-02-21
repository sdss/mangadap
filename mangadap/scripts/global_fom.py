#!/usr/bin/env python3
import os
import time
import glob
import numpy

from matplotlib import pyplot, colors, ticker, cm

from scipy.special import erf
from scipy import interpolate

from astropy.io import fits
from astropy.stats import median_absolute_deviation

from mangadap.util.fileio import init_record_array, rec_to_fits_type, rec_to_fits_col_dim

#-----------------------------------------------------------------------------

def _dtype():
    return [ ('X', numpy.float64),
             ('L', numpy.float64),
             ('M', numpy.float64),
             ('H', numpy.float64) ]

def twod_bin(x, y, x0, dx, nx, y0, dy, ny):
    bini = ((x-x0)/dx).astype(int)
    bini[ (bini < 0) | (bini >= nx) ] = -1
    binj = ((y-y0)/dy).astype(int)
    binj[ (binj < 0) | (binj >= ny) ] = -1

    binij = numpy.ma.MaskedArray(bini*ny + binj, mask=(bini < 0) | (binj < 0))

    uniq, cnt = numpy.unique(binij.compressed(), return_counts=True)
    n = numpy.zeros((nx,ny),dtype=int)
    i = uniq//ny
    j = uniq - i*ny
    n[i,j] = cnt

    return n

def medstats(x):
    srt = numpy.sort(x)
    samp = (numpy.array([0.05, 0.5, 0.95])*len(x)).astype(int)
    return srt[samp]

def medstats_bin(x, y, x0, dx, nx):
    bini = ((x-x0)/dx).astype(int)
    bini[ (bini < 0) | (bini >= nx) ] = -1

    db = init_record_array(nx, _dtype())

    for i in range(nx):
        print('{0}/{1}'.format(i+1, nx), end='\r')
        indx = bini == i
        if not numpy.any(indx):
            continue
        db['X'][i] = numpy.median(x[indx])
        db['L'][i], db['M'][i], db['H'][i] = medstats(y[indx])
    return db

def get_data_sc(files, names, pltifu=None):
    n = len(names)
    d = [numpy.empty(0, dtype=float)]
    for i in range(1,n):
        d += [numpy.empty(0, dtype=float)]
    plate = None if pltifu is None else int(pltifu.split('-')[0])
    for f in files:
        if plate is not None and str(plate) not in f:
            continue
        hdu = fits.open(f)
        plateifu = numpy.unique(hdu[1].data['PLATEIFU']) if plate is None else [pltifu]
        for pi in plateifu:
            pi_indx = hdu[1].data['PLATEIFU'] == pi
            uniq, indx = numpy.unique(hdu[1].data['BINID'][pi_indx,0], return_index=True)
            if uniq[0] == -1:   # This shouldn't happen, but just in case ...
                indx = indx[1:]
            for i in range(n):
                if names[i] == 'snr':
                    d[i] = numpy.append(d[i], hdu[1].data['SNR'][pi_indx,1][indx])
                elif names[i] == 'chi':
                    d[i] = numpy.append(d[i], hdu[1].data['SFM'][pi_indx,3][indx])
                elif names[i] == 'frm':
                    d[i] = numpy.append(d[i], hdu[1].data['SFM'][pi_indx,1][indx])
                elif names[i] == 'g68':
                    d[i] = numpy.append(d[i], hdu[1].data['SFMG'][pi_indx,2,0][indx])
                elif names[i] == 'g99':
                    d[i] = numpy.append(d[i], hdu[1].data['SFMG'][pi_indx,2,2][indx])
        hdu.close()
        del hdu
    return d

def get_data_el(files, names, pltifu=None):
    n = len(names)
    d = [numpy.empty(0, dtype=float)]
    for i in range(1,n):
        d += [numpy.empty(0, dtype=float)]
    plate = None if pltifu is None else int(pltifu.split('-')[0])
    for f in files:
        if plate is not None and str(plate) not in f:
            continue
        hdu = fits.open(f)
        indx = numpy.ones(hdu[1].data['SNR'].shape[0], dtype=bool) if plate is None else \
                    hdu[1].data['PLATEIFU'] == pltifu
        for i in range(n):
            if names[i] == 'snr':
                d[i] = numpy.append(d[i], hdu[1].data['SNR'][indx,0])
            elif names[i] == 'chi':
                d[i] = numpy.append(d[i], hdu[1].data['EFM'][indx,3])
            elif names[i] == 'frm':
                d[i] = numpy.append(d[i], hdu[1].data['EFM'][indx,1])
            elif names[i] == 'g68':
                d[i] = numpy.append(d[i], hdu[1].data['EFMG'][indx,2,0])
            elif names[i] == 'g99':
                d[i] = numpy.append(d[i], hdu[1].data['EFMG'][indx,2,2])
        hdu.close()
        del hdu
    return d

def init_ax(fig, pos, facecolor=None):
    ax = fig.add_axes(pos, facecolor=facecolor)
    ax.minorticks_on()
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
#    ax.tick_params(which='both', direction='in') #, top='on', right='on')
    ax.grid(True, which='major', color='0.8', zorder=1, linestyle='-', lw=0.5)
    return ax


def make_log_grid(xbins, ybins):
    dx = numpy.mean(numpy.diff(numpy.log10(xbins)))
    dy = numpy.mean(numpy.diff(numpy.log10(ybins)))
    X, Y = numpy.meshgrid(xbins[:-1]+dx*xbins[:-1]*numpy.log(10)/2,
                          ybins[:-1]+dy*ybins[:-1]*numpy.log(10)/2)
    return dx, dy, X, Y


def fill_log_grid(xname, yname, xbins, ybins, ofile, sc, overwrite, snlim=None):
    if os.path.isfile(ofile) and not overwrite:
        return

    files = glob.glob('chi2_tests/*.fits')

    dx, dy, X, Y = make_log_grid(xbins, ybins)

#    pltifu = '9193-12704'
#    snr, chi, g68, g99 = get_data_sc(files, pltifu=pltifu) \
#                            if sc else get_data_el(files, pltifu=pltifu)
##    pyplot.scatter(chi, g68/m[0], marker='.', s=20, lw=0, color='C3')
#    pyplot.scatter(chi, numpy.ma.divide(g99,g68) * m[0]/m[2], marker='.', s=20, lw=0, color='C3')
##    pyplot.scatter(snr, chi, marker='.', s=20, lw=0, color='C3')
#    pyplot.xscale('log')
#    pyplot.yscale('log')
#    pyplot.show()
#    exit()

    print('getting data')
    x, y = get_data_sc(files, [xname, yname]) if sc else get_data_el(files, [xname, yname])

    print('binning data')
    x0 = numpy.log10(xbins[0])
    y0 = numpy.log10(ybins[0])

    _x = numpy.ma.log10(x)
    _y = numpy.ma.log10(y)

    mask = _x.mask | _y.mask
    _x[mask] = numpy.ma.masked
    _y[mask] = numpy.ma.masked

    db = medstats_bin(_x.compressed(), _y.compressed(), x0, dx*2, nbins//2)
    binned = twod_bin(_x.compressed(), _y.compressed(), x0, dx, nbins, y0, dy, nbins)

    binned = twod_bin(_x.compressed(), _y.compressed(), x0, dx, nbins, y0, dy, nbins)
    print('writing: {0}'.format(ofile))
    hdr = fits.Header()
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.ImageHDU(data=binned),
                  fits.BinTableHDU.from_columns([
                        fits.Column(name=n, format=rec_to_fits_type(db[n]),
                                    dim=rec_to_fits_col_dim(db[n]), array=db[n])
                                        for n in db.dtype.names], name='GROW')
                 ]).writeto(ofile, overwrite=overwrite)

#    binned = numpy.expand_dims(twod_bin(numpy.log10(x), numpy.log10(y), x0, dx, nbins, y0, dy,
#                               nbins), axis=2)
#    for i in range(len(sn)-1):
#        indx = (snr > sn[i]) & (snr < sn[i+1])
#        hdr['SN{0}'.format(i+1)] = ('{0},{1}'.format(sn[i],sn[i+1]), 'SN range {0}'.format(i+1))
#        binned = numpy.append(binned, numpy.expand_dims(twod_bin(numpy.log10(chi[indx]),
#                                                                 numpy.log10(d[indx]), x0, dx,
#                                                                 nbins, y0, dy, nbins),
#                                                        axis=2), axis=2)
#        print('SN{0} binned'.format(i+1))

def global_fom_data(analysis_path, daptype, ofile):
    pass
    
    

#-----------------------------------------------------------------------------

def junk():
    t = time.perf_counter()

    x = numpy.linspace(0,5,100)
    g = (erf(x/numpy.sqrt(2)) - erf(-x/numpy.sqrt(2)))/2.
    interp = interpolate.interp1d(g, x)
    m = interp([0.68, 0.95, 0.99])
    print(m[2]/m[0])

    nbins = 256
    chi_lim = [0.5,  1000]
    grw_lim = [0.4,     7]
    snr_lim = [0.9,   200]
    frm_lim = [0.008, 200]

    sn = [ 1, 5, 10, 20, 40, 80 ]

    snrbins = numpy.logspace(*numpy.log10(snr_lim),nbins+1)
    chibins = numpy.logspace(*numpy.log10(chi_lim),nbins+1)
    grwbins = numpy.logspace(*numpy.log10(grw_lim),nbins+1)
    frmbins = numpy.logspace(*numpy.log10(frm_lim),nbins+1)

    sc = False #True
    overwrite = False #True

    ofile = 'snr_frm_sc.fits' if sc else 'snr_frm_el.fits'
    fill_log_grid('snr', 'frm', snrbins, frmbins, ofile, sc, overwrite)

    hdu = fits.open(ofile)
    binned = hdu[1].data
    db = hdu[2].data
    hdu.close()

    print('plotting data')
    n_lim = [1, numpy.amax(binned)] #10000]
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.1, 0.1, 0.8, 0.8], facecolor='0.95')
    ax.set_xlim(snr_lim)
    ax.set_ylim(frm_lim)
    ax.set_xscale('log')
    ax.set_yscale('log')
#    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
#    ax.yaxis.set_major_locator(ticker.LogLocator(base=10,subs=(1,2,4,8,)))
    ax.text(0.5, -0.07, r'S/N$_g$', horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.text(-0.1, 0.5, r'fRMS', horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, rotation='vertical')

#    cmap = cm.get_cmap('inferno')
#    cnorm = colors.LogNorm(vmin=sn[1]/1.1, vmax=sn[-1]*1.1)
#    print(sn)

    tot_binned = numpy.ma.MaskedArray(binned[:,:], mask=numpy.invert(binned[:,:]>0)).T
#    tot_binned = numpy.ma.MaskedArray(binned[:,:,0], mask=numpy.invert(binned[:,:,0]>0)).T
    print(numpy.amax(tot_binned))
    img = ax.pcolormesh(snrbins, frmbins, tot_binned,
                        norm=colors.LogNorm(vmin=n_lim[0], vmax=n_lim[1]),
                        cmap='viridis', zorder=4, lw=0, rasterized=True)

    x = numpy.power(10, db['X'])
    l = numpy.power(10, db['L'])
    m = numpy.power(10, db['M'])
    h = numpy.power(10, db['H'])

    indx = (x > 3) & (x < 80)

    ax.plot(x[indx], l[indx], lw=2.0, color='k', linestyle='--', zorder=7, alpha=0.7)
    ax.plot(x[indx], h[indx], lw=2.0, color='k', linestyle='--', zorder=7, alpha=0.7)
    ax.plot(x[indx], m[indx], lw=2.0, color='k', linestyle='-', zorder=7, alpha=0.7)

    pyplot.show()
    exit()

    nsrt = numpy.sort(tot_binned.compressed())[::-1]
    cumulative = numpy.cumsum(nsrt)
    interp = interpolate.interp1d(cumulative, nsrt)
    lev = interp([0.99*cumulative[-1]])
    print(lev)
#    ax.contour(X, Y, tot_binned, levels=lev, linewidths=1, zorder=5, colors='0.5')
    for i in range(1,len(sn)):
        tot_binned = numpy.ma.MaskedArray(binned[:,:,i], mask=numpy.invert(binned[:,:,i]>0)).T
        nsrt = numpy.sort(tot_binned.compressed())[::-1]
        cumulative = numpy.cumsum(nsrt)
        interp = interpolate.interp1d(cumulative, nsrt)
        lev = interp([0.99*cumulative[-1]]) if i < len(sn)-1 else interp([0.95*cumulative[-1]])
        print(lev)
        print(sn[i])
        print(cmap(cnorm(sn[i])))
        ax.contour(X, Y, tot_binned, levels=lev, linewidths=1, zorder=5,
                   colors=[cmap(cnorm(sn[i]))])
#                   cmap=cmap, norm=colors.LogNorm(vmin=md_lim[0]/3, vmax=md_lim[1]/3))

    pltifu = '9193-12704'
    snr, chi, g68, g99 = get_data_sc(files, pltifu=pltifu) \
                            if sc else get_data_el(files, pltifu=pltifu)
    d = g99/g68 * m[0]/m[2]
    for i in range(len(sn)-1):
        indx = (snr > sn[i]) & (snr < sn[i+1])
        ax.scatter(chi[indx], d[indx], marker='.', s=40, lw=0.5, edgecolor='k', color=cmap(cnorm(sn[i+1])), zorder=6)
    pyplot.show()
    exit()

    ngood = 0
    ntot = 0
    for f in files[:50]:

        hdu = fits.open(f)
#        d = hdu[1].data['EFMG'][:,2,2]/m[2] - hdu[1].data['EFMG'][:,2,1]/m[1]

        indx = hdu[1].data['SNR'][:,0] > 40
        d = hdu[1].data['EFMG'][indx,2,2]/hdu[1].data['EFMG'][indx,2,0] * m[0]/m[2]

#        good = numpy.sum((d > -0.15) & (d < 1) & (hdu[1].data['EFM'][:,3] > 3))
#        print(f)
#        print('{0:.3f}'.format(good/hdu[1].data['EFMG'].shape[0]))
#        ntot += hdu[1].data['EFMG'].shape[0]
#        ngood += good

        pyplot.scatter(hdu[1].data['EFM'][indx,3], d, marker='.', lw=0, s=20, alpha=0.5)
#        pyplot.scatter(hdu[1].data['EFM'][:,3], d, marker='.', lw=0, s=20, alpha=0.5)
#    pyplot.scatter(hdu[1].data['EFM'][:,3], hdu[1].data['EFMG'][:,2,2]/m[2], marker='.', lw=0, s=20,
#                   alpha=0.5)
#    pyplot.scatter(hdu[1].data['EFM'][:,3], hdu[1].data['EFMG'][:,2,1]/m[1], marker='.', lw=0, s=20,
#                   alpha=0.5)
#    pyplot.scatter(hdu[1].data['EFM'][:,3], hdu[1].data['EFMG'][:,2,0]/m[0], marker='.', lw=0, s=20,
#                   alpha=0.5)
#    pyplot.plot([0.5,30], [0.5,30], color='C3')
#    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel(r'$\chi^2_\nu$')
#    pyplot.ylabel(r'Normalized Growth')
#    print('{0:.3f}'.format(ngood/ntot))
    pyplot.show()
    exit()

    pyplot.scatter(hdu[1].data['SNR'][:,0], hdu[1].data['EFMG'][:,2,2]/m[2], marker='.', lw=0, s=20,
                   alpha=0.5)
    pyplot.scatter(hdu[1].data['SNR'][:,0], hdu[1].data['EFMG'][:,2,1]/m[1], marker='.', lw=0, s=20,
                   alpha=0.5)
    pyplot.scatter(hdu[1].data['SNR'][:,0], hdu[1].data['EFMG'][:,2,0]/m[0], marker='.', lw=0, s=20,
                   alpha=0.5)
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel(r'S/N$_g$')
    pyplot.ylabel(r'Normalized Growth')
    pyplot.show()

    pyplot.scatter(hdu[1].data['SNR'][:,1], hdu[1].data['SFM'][:,3], marker='.', lw=0, s=20,
                   alpha=0.5)
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel(r'S/N$_g$')
    pyplot.ylabel(r'Stellar $\chi^2_\nu$')
    pyplot.show()

    pyplot.scatter(hdu[1].data['SNR'][:,0], hdu[1].data['EFM'][:,3], marker='.', lw=0, s=20,
                   alpha=0.5)
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel(r'S/N$_g$')
    pyplot.ylabel(r'Full $\chi^2_\nu$')
    pyplot.show()

    pyplot.scatter(hdu[1].data['EFM'][:,3], hdu[1].data['EFMG'][:,2,2]/m[2], marker='.', lw=0, s=20,
                   alpha=0.5)
#    pyplot.scatter(hdu[1].data['EFM'][:,3], hdu[1].data['EFMG'][:,2,1]/m[1], marker='.', lw=0, s=20,
#                   alpha=0.5)
#    pyplot.scatter(hdu[1].data['EFM'][:,3], hdu[1].data['EFMG'][:,2,0]/m[0], marker='.', lw=0, s=20,
#                   alpha=0.5)
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel(r'$\chi^2_\nu$')
    pyplot.ylabel(r'Normalized Growth')
    pyplot.show()
    exit()




    pyplot.plot(x,1-g)
    pyplot.yscale('log')
    pyplot.show()

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



    pyplot.scatter(hdu[1].data['SNR'][:,1], hdu[1].data['SFMG'][:,2,0], marker='.', s=20, lw=0, alpha=0.5)
    pyplot.scatter(hdu[1].data['SNR'][:,1], hdu[1].data['SFMG'][:,2,1], marker='.', s=20, lw=0, alpha=0.5)
    pyplot.scatter(hdu[1].data['SNR'][:,1], hdu[1].data['SFMG'][:,2,2], marker='.', s=20, lw=0, alpha=0.5)
    pyplot.plot([0.1,100], [2.57583409,2.57583409], color='C2', linestyle=':')
    pyplot.plot([0.1,100], [1.96034685,1.96034685], color='C1', linestyle=':')
    pyplot.plot([0.1,100], [0.99472642,0.99472642], color='C0', linestyle=':')
    pyplot.xscale('log')
    pyplot.show()


if __name__ == '__main__':
    t = time.perf_counter()

    parser = argparse.ArgumentParser()

    parser.add_argument('--drpver', type=str, help='DRP version', default=None)
    parser.add_argument('--dapver', type=str, help='DAP version', default=None)
    parser.add_argument('--dap_src', type=str, help='Top-level directory with the DAP source code;'
                        ' defaults to $MANGADAP_DIR', default=None)
    parser.add_argument("--redux_path", type=str, help="main DRP output path", default=None)
    parser.add_argument("--analysis_path", type=str, help="main DAP output path", default=None)

    parser.add_argument("--plan_file", type=str, help="parameter file with the MaNGA DAP "
                        "execution plan to use instead of the default" , default=None)

    parser.add_argument('--daptype', type=str, help='DAP processing type', default=None)
    parser.add_argument('--normal_backend', dest='bgagg', action='store_false', default=True)

    parser.add_argument('--force', type=bool, action='store_true', default=False,
                        help='Force the reconstruction of the global FOM data file.')

    arg = parser.parse_args()

    if arg.bgagg:
        pyplot.switch_backend('agg')

    # Set the paths
    redux_path = defaults.default_redux_path(drpver=arg.drpver) \
                        if arg.redux_path is None else arg.redux_path
    analysis_path = defaults.default_analysis_path(drpver=arg.drpver, dapver=arg.dapver) \
                            if arg.analysis_path is None else arg.analysis_path

    daptypes = []
    if arg.daptype is None:
        plan_file = defaults.default_dap_plan_file(drpver=arg.drpver, dapver=arg.dapver,
                                          analysis_path=arg.analysis_path) \
                                            if arg.plan_file is None else arg.plan_file
        analysisplan = AnalysisPlanSet.from_par_file(plan_file)
#        daptypes = [ defaults.default_dap_method(plan=p) for p in analysisplan ]
        for p in analysisplan:
            bin_method = SpatiallyBinnedSpectra.define_method(p['bin_key'])
            sc_method = StellarContinuumModel.define_method(p['continuum_key'])
            el_method = EmissionLineModel.define_method(p['elfit_key'])
            daptypes += [defaults.default_dap_method(bin_method['key'],
                                                     sc_method['fitpar']['template_library_key'],
                                                     el_method['continuum_tpl_key'])]
    else:
        daptypes = [arg.daptype]

    for daptype in daptypes:
        oroot = os.path.join(analysis_path, daptype, 'qa')
        ofile = os.path.join(oroot, 'global_fom_{0}.fits.gz'.format(daptype))
        if not os.path.isfile(ofile) or arg.force:
            global_fom_data(analysis_path, daptype, ofile)

        oplot = os.path.join(oroot, 'global_fom_{0}.pdf'.format(daptype))
        global_fom_plots(analysis_path, daptype)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


