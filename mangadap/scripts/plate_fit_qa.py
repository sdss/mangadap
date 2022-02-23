import os
import time
import argparse
import warnings
import numpy

from matplotlib import pyplot, ticker, rc

from astropy.io import fits

from mangadap.dapfits import DAPQualityBitMask
from mangadap.config import defaults
from mangadap.proc.util import growth_lim
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.util.fileio import channel_dictionary

from mangadap.scripts import scriptbase


def init_ax(fig, pos, facecolor=None, grid=False):
    ax = fig.add_axes(pos, facecolor=facecolor)
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    if grid:
        ax.grid(True, which='major', color='0.8', zorder=1, linestyle='-', lw=0.5)
    return ax


def is_critical(dapver, analysis_path, daptype, plt):

    ifus = [1901, 1902, 3701, 3702, 3703, 3704, 6101, 6102, 6103, 6104, 9101, 9102,
            12701, 12702, 12703, 12704, 12705]

    dapqualbm = DAPQualityBitMask()

    critical = numpy.zeros(len(ifus), dtype=bool)

    for i,ifu in enumerate(ifus):
        apath = os.path.join(analysis_path, daptype, str(plt), str(ifu))
        maps_file = os.path.join(apath, 'manga-{0}-{1}-MAPS-{2}.fits.gz'.format(plt,ifu,daptype))

        if not os.path.isfile(maps_file):
            continue

        hdu = fits.open(maps_file)
        critical[i] = dapqualbm.flagged(hdu[0].header['dapqual'], 'CRITICAL')

    return ifus, critical


def compile_data(dapver, analysis_path, daptype, plt):

    ifus = [1901, 1902, 3701, 3702, 3703, 3704, 6101, 6102, 6103, 6104, 9101, 9102,
            12701, 12702, 12703, 12704, 12705]

    sc_pltifu = []
    sc_snrg = []
    sc_frms = []
    sc_rchi = []

    el_pltifu = []
    el_snrg = []
    el_frms = []
    el_rchi = []

    for ifu in ifus:
        apath = os.path.join(analysis_path, daptype, str(plt), str(ifu))
        maps_file = os.path.join(apath, 'manga-{0}-{1}-MAPS-{2}.fits.gz'.format(plt,ifu,daptype))

        if not os.path.isfile(maps_file):
            continue

        hdu = fits.open(maps_file)

        # Get the stellar data
#        uniq, indx = map(lambda x: x[1:], numpy.unique(hdu['BINID'].data[1,:,:].ravel(),
#                                                       return_index=True))
        uniq, indx = numpy.unique(hdu['BINID'].data[1,:,:].ravel(), return_index=True)
        if uniq[0] == -1:
            uniq = uniq[1:]
            indx = indx[1:]

        mask = hdu['STELLAR_VEL_MASK'].data.ravel()[indx] > 0
        indx = indx[numpy.invert(mask)]

        if len(indx) > 0:
            sc_pltifu += [ '{0}-{1}'.format(plt, ifu) ]
            sc_snrg += [ hdu['BIN_SNR'].data.ravel()[indx] ]
            sc_frms += [ hdu['STELLAR_FOM'].data[1,:,:].ravel()[indx] ]
            sc_rchi += [ hdu['STELLAR_FOM'].data[2,:,:].ravel()[indx] ]

        # Get the emission-line data
        eml = channel_dictionary(hdu, 'EMLINE_GFLUX')
#        uniq, indx = map(lambda x: x[1:], numpy.unique(hdu['BINID'].data[3,:,:].ravel(),
#                                                       return_index=True))
        uniq, indx = numpy.unique(hdu['BINID'].data[3,:,:].ravel(), return_index=True)
        if uniq[0] == -1:
            uniq = uniq[1:]
            indx = indx[1:]
        mask = hdu['EMLINE_GFLUX_MASK'].data[eml['Ha-6564'],:,:].ravel()[indx] > 0
        indx = indx[numpy.invert(mask)]

        if len(indx) > 0:
            el_pltifu += [ '{0}-{1}'.format(plt, ifu) ]
            el_snrg += [ hdu['SPX_SNR'].data.ravel()[indx] ]
            el_frms += [ hdu['EMLINE_FOM'].data[1,:,:].ravel()[indx] ]
            el_rchi += [ hdu['EMLINE_FOM'].data[2,:,:].ravel()[indx] ]

    return ifus, numpy.array(sc_pltifu), sc_snrg, sc_frms, sc_rchi, numpy.array(el_pltifu), \
                el_snrg, el_frms, el_rchi


def plate_fit_qa(dapver, analysis_path, daptype, plt):

    # Make sure the output directory exists
    plate_qa_dir = os.path.join(analysis_path, daptype, str(plt), 'qa')
    if not os.path.isdir(plate_qa_dir):
        os.makedirs(plate_qa_dir)

    ofile = os.path.join(plate_qa_dir, '{0}-fitqa.png'.format(plt))
#    ofile = None

    # Determine if the analysis is a CRITICAL failure
    c_ifus, critical = is_critical(dapver, analysis_path, daptype, plt)
    pltifu = ['{0}-{1}'.format(plt, ifu) for ifu in c_ifus]

    # Grab the data
    d_ifus, sc_pltifu, sc_snrg, sc_frms, sc_rchi, el_pltifu, el_snrg, el_frms, el_rchi \
            = compile_data(dapver, analysis_path, daptype, plt)

    # Check if there's anything to plot
    if len(sc_pltifu) == 0 & len(el_pltifu) == 0:
        warnings.warn('No PLATEIFU analysis complete for plate={0}, daptype={1}.'.format(
                                                                                    plt, daptype))
        return

    # Basic coding check
    assert numpy.array_equal(c_ifus, d_ifus), 'Coding error'

    # Make the plot
    rc('font', size=10)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    axsc = init_ax(fig, [0.08, 0.54, 0.45, 0.45])
    axec = init_ax(fig, [0.54, 0.54, 0.45, 0.45])
    axsf = init_ax(fig, [0.08, 0.08, 0.45, 0.45])
    axef = init_ax(fig, [0.54, 0.08, 0.45, 0.45])

#    chi = numpy.concatenate(sc_rchi+el_rchi)
#    chi_lim = numpy.power(10, growth_lim(numpy.ma.log10(chi).compressed(), 0.999, fac=1.5))
#
#    frms = numpy.append(sc_frms+el_frms)
#    frms_lim = numpy.power(10, growth_lim(numpy.ma.log10(frms).compressed(), 0.999, fac=1.5))
#
#    snr = numpy.append(sc_snrg+el_snrg)
#    snr_lim = numpy.power(10, growth_lim(numpy.ma.log10(snr).compressed(), 0.999, fac=1.5))

    chi_lim = [0.3, 100]
    frms_lim = [0.005, 500]
    snr_lim = [0.8, 300]

    axsc.set_xlim(snr_lim)
    axec.set_xlim(snr_lim)
    axsf.set_xlim(snr_lim)
    axef.set_xlim(snr_lim)

    axsc.set_ylim(chi_lim)
    axec.set_ylim(chi_lim)

    axsf.set_ylim(frms_lim)
    axef.set_ylim(frms_lim)

    axsc.set_xscale('log')
    axsc.set_yscale('log')
    axec.set_xscale('log')
    axec.set_yscale('log')
    axsf.set_xscale('log')
    axsf.set_yscale('log')
    axef.set_xscale('log')
    axef.set_yscale('log')

    axsc.tick_params(axis='x', which='both', direction='in')
    axec.tick_params(axis='x', which='both', direction='in')
    axec.tick_params(axis='y', which='both', direction='in')
    axef.tick_params(axis='y', which='both', direction='in')

    axsc.xaxis.set_major_formatter(ticker.NullFormatter())
    axec.xaxis.set_major_formatter(ticker.NullFormatter())
    axec.yaxis.set_major_formatter(ticker.NullFormatter())
    axef.yaxis.set_major_formatter(ticker.NullFormatter())

    # Plot the data
    plotted = numpy.zeros(len(pltifu), dtype=bool)
    for i in range(len(pltifu)):
        # Set the point type based on the critical flag and the number
        # of pltifus plotted
        s = 20 if critical[i] else 40
        c = 'C{0}'.format(i % 10)
        lw = 0 if i < 10 else 0.5

        # Find the correct data
        sc_indx = [False] if len(sc_pltifu) == 0 else numpy.array(sc_pltifu) == pltifu[i]
        el_indx = [False] if len(el_pltifu) == 0 else numpy.array(el_pltifu) == pltifu[i]

        if numpy.sum(sc_indx) == 0 and numpy.sum(el_indx) == 0:
            continue

        plotted[i] = True

        # Plot the stellar data
        if numpy.sum(sc_indx) == 1:
            j = numpy.where(sc_indx)[0][0]
            axsc.scatter(sc_snrg[j], sc_rchi[j], marker='.', edgecolor='k', color=c, s=s, lw=lw,
                         alpha=0.7, zorder=3)
            axsf.scatter(sc_snrg[j], sc_frms[j], marker='.', edgecolor='k', color=c, s=s, lw=lw,
                         alpha=0.7, zorder=3)
        if numpy.sum(el_indx) == 1:
            j = numpy.where(sc_indx)[0][0]
            axec.scatter(el_snrg[j], el_rchi[j], marker='.', edgecolor='k', color=c, s=s, lw=lw,
                         alpha=0.7, zorder=3, label=pltifu[i])
            axef.scatter(el_snrg[j], el_frms[j], marker='.', edgecolor='k', color=c, s=s, lw=lw,
                         alpha=0.7, zorder=3)
        else:
            axec.scatter([], [], marker='.', edgecolor='k', color=c, s=s, lw=lw,
                         alpha=0.7, zorder=3, label=pltifu[i])
        
    axsc.text(-0.13, 0.5, r'$\chi^2_\nu$', ha='center', va='center', transform=axsc.transAxes,
              rotation='vertical', fontsize=12)
    axsf.text(-0.13, 0.5, r'Fractional RMS', ha='center', va='center', transform=axsf.transAxes,
              rotation='vertical', fontsize=12)
    axsf.text(0.5, -0.13, r'S/N$_g$', ha='center', va='center', transform=axsf.transAxes,
              fontsize=12)
    axef.text(0.5, -0.13, r'S/N$_g$', ha='center', va='center', transform=axef.transAxes,
              fontsize=12)

    l = axec.legend(fontsize='x-small', ncol=2, loc=2)
    if numpy.any(critical[plotted]):
        indx = numpy.where(critical[plotted])[0]
        t = l.get_texts()
        for i in indx:
            t[i].set_color('r')

    if ofile is None:
        pyplot.show()
    else:
        print('Writing: {0}'.format(ofile))
        fig.canvas.print_figure(ofile, bbox_inches='tight', dpi=100)
    fig.clear()
    pyplot.close(fig)


class PlateFitQA(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Construct per-plate QA plots', width=width)

        parser.add_argument('plate', type=int, help='plate ID to process')

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

        for daptype in daptypes:
            plate_fit_qa(args.dapver, analysis_path, daptype, args.plate)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

