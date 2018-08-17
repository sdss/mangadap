#!/usr/bin/env python3

import os
import time
import numpy

from matplotlib import pyplot, ticker, rc

from astropy.io import fits

from mangadap.proc.util import growth_lim
from mangadap.util.bitmask import BitMask

#-----------------------------------------------------------------------------
def init_ax(fig, pos):
    ax = fig.add_axes(pos) #, facecolor='0.9')
    ax.minorticks_on()
#    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
    ax.minorticks_on()
    ax.tick_params(which='major', length=6) #, direction='in')
    ax.tick_params(which='minor', length=3) #, direction='in')
    return ax


def radial_coverage_histogram(dapall, plot_file=None):

    # Get the range
    arcsec_rng = growth_lim(dapall['RCOV90'], 0.99, fac=1.1)
    rore = numpy.ma.divide(dapall['RCOV90'], dapall['NSA_ELPETRO_TH50_R'])
    rore[dapall['NSA_ELPETRO_TH50_R'] < 0] = numpy.ma.masked
    rore_rng = growth_lim(rore, 0.99, fac=1.1)

    # Determine the sample
#    sdssMaskbits = os.path.join(os.environ['IDLUTILS_DIR'], 'data', 'sdss', 'sdssMaskbits.par')
    sdssMaskbits = os.path.join(os.environ['MANGADAP_DIR'], 'data', 'sdss', 'sdssMaskbits.par')
    targ1bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET1')

    ancillary = dapall['MNGTARG3'] > 0
    primaryplus = targ1bm.flagged(dapall['MNGTARG1'],
                                  flag=['PRIMARY_v1_2_0', 'COLOR_ENHANCED_v1_2_0'])
    secondary = targ1bm.flagged(dapall['MNGTARG1'], flag='SECONDARY_v1_2_0')
    other = numpy.invert(ancillary) & numpy.invert(primaryplus) & numpy.invert(secondary)

    indx = dapall['RCOV90'] < 1
    print('Num. of observations with 90% coverage <1 arcsec: {0}'.format(numpy.sum(indx)))
    print(dapall['PLATEIFU'][indx])

    # Instantiate the figure
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    nlim = [0.8, 1000]

    # Show RCOV90 in terms of arcseconds for each IFU
    ax = init_ax(fig, [0.1, 0.58, 0.8, 0.4])
    ax.set_ylim(nlim)
    ax.set_yscale('log', nonposy='clip')
    ifu = [ 19, 37, 61, 91, 127 ]
    for i,f in enumerate(ifu):
        indx = (dapall['IFUDESIGN']/100).astype(int) == f
        ax.hist(dapall['RCOV90'][indx], bins=50, range=arcsec_rng, alpha=0.5,
                histtype='stepfilled', color='C{0}'.format(i), zorder=2)
        ax.hist(dapall['RCOV90'][indx], bins=50, range=arcsec_rng,
                histtype='step', color='C{0}'.format(i), zorder=3)
        medp = numpy.median(dapall['RCOV90'][indx])
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
    ax.set_yscale('log', nonposy='clip')
    ax.hist(rore[primaryplus], bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
            color='C0', zorder=2)
    ax.hist(rore[primaryplus], bins=50, range=rore_rng, histtype='step',
            color='C0', zorder=3)
    medp = numpy.median(rore[primaryplus])
    ax.plot([medp, medp], nlim, color='C0', zorder=4, linestyle='--')
    ax.hist(rore[secondary], bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
            color='C1', zorder=2)
    ax.hist(rore[secondary], bins=50, range=rore_rng, histtype='step',
            color='C1', zorder=3)
    medp = numpy.median(rore[secondary])
    ax.plot([medp, medp], nlim, color='C1', zorder=4, linestyle='--')
    ax.hist(rore[other], bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
            color='C2', zorder=2)
    ax.hist(rore[other], bins=50, range=rore_rng, histtype='step',
            color='C2', zorder=3)
    ax.hist(rore[ancillary], bins=50, range=rore_rng, alpha=0.5, histtype='stepfilled',
            color='C3', zorder=2)
    ax.hist(rore[ancillary], bins=50, range=rore_rng, histtype='step',
            color='C3', zorder=3)
    ax.text(-0.1, 0.5, 'N', ha='center', va='center', transform=ax.transAxes, rotation='vertical')
    ax.text(0.5, -0.12, r'$R_{90}/R_{eff}$', ha='center', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.95, 'Primary+', color='C0', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.90, 'Secondary', color='C1', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.85, 'Other', color='C2', ha='left', va='center', transform=ax.transAxes)
    ax.text(0.05, 0.80, 'Ancillary', color='C3', ha='left', va='center', transform=ax.transAxes)

    if plot_file is None:
        pyplot.show()
    else:
        fig.canvas.print_figure(plot_file, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    dapall_file = '/Users/westfall/Work/MaNGA/testing/dr15/dapall-v2_4_3-2.2.1.fits'
    dapall_hdu = fits.open(dapall_file)

    radial_coverage_histogram(dapall_hdu['DAPALL'].data,
                              plot_file='../ms/figs/radial_coverage.pdf')

    print('Elapsed time: {0} seconds'.format(time.clock() - t))


