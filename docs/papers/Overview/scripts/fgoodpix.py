#!/usr/bin/env python3

import os
import time
import numpy
import glob

from astropy.io import fits
from matplotlib import pyplot, rc, cm, colors, ticker

os.environ['MANGADRP_VER'] = 'v2_4_3'
os.environ['MANGA_SPECTRO_REDUX'] = '/Volumes/MaNGA/redux'
os.environ['MANGADAP_VER'] = '2.2.1'
os.environ['MANGA_SPECTRO_ANALYSIS'] = '/Volumes/MaNGA/analysis'

from mangadap.drpfits import DRPFits

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()

    font = { 'size' : 16 }
    rc('font', **font)

    cmap = cm.get_cmap('magma')

    ifulist = [  1901,  1902,  3701,  3702,  3703, 3704, 6101, 6102, 6103, 6104, 9101, 9102,
                12701, 12702, 12703, 12704, 12705 ]

#    c = [ 'C0','C0','C1','C1','C1','C1','C2','C2','C2','C2','C3','C3','C4','C4','C4','C4','C4']

    normed_color = colors.Normalize(vmin=0.5, vmax=6.0)
    c = [ cmap(normed_color(1)) ]*2 + [ cmap(normed_color(2)) ]*4 + [ cmap(normed_color(3)) ]*4 \
            + [ cmap(normed_color(4)) ]*2 + [ cmap(normed_color(5)) ]*5

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = fig.add_axes([0.12,0.1,0.8,0.6])
    ax.minorticks_on()
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
#    ax.tick_params(which='both', top=True, right=True)
    ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-', lw=0.5)

#    ylim = [0.11, 0.59]
    ylim = [0.5, 1.02]
    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim(ylim)

#    for ifu,clr in zip([ifulist[1],ifulist[-2]],[c[1],c[-2]]):
    for ifu,clr in zip(ifulist,c):
        drpf = DRPFits(7495, ifu, 'CUBE', drpver='v2_4_3', read=True)
        print(drpf.file_path())
        flux = drpf.copy_to_masked_array(flag=['DONOTUSE', 'FORESTAR'])
        fgoodpix = numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux)),axis=1) / flux.shape[1]

        if ifu == 12704:
            fgoodpix_saved = fgoodpix.copy().reshape(drpf.spatial_shape)

#        print(fgoodpix.shape)
        srt = numpy.argsort(fgoodpix)
        _fgoodpix = numpy.ma.MaskedArray(fgoodpix[srt], mask=numpy.invert(fgoodpix[srt] > 0))
        omega = 1-(numpy.arange(len(srt), dtype=float)+1)/len(srt)

        omega /= numpy.amax(omega[numpy.invert(_fgoodpix.mask)])

#        print(fgoodpix[srt])
#        ax.step(_fgoodpix, 1-(numpy.arange(len(srt), dtype=float)+1)/len(srt),
#                    color=clr, lw=1.5, zorder=3)
        ax.step(_fgoodpix, omega, color=clr, lw=1.5, zorder=3)

    ax.text(0.50, -0.1, r'$\delta\Lambda$', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.12, 0.50, r'$\delta\Omega (>\delta\Lambda)$', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')

    ax.text(0.05, 0.32, r'7495-19xx', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes, color=cmap(normed_color(1)))
    ax.text(0.05, 0.26, r'7495-37xx', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes, color=cmap(normed_color(2)))
    ax.text(0.05, 0.20, r'7495-61xx', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes, color=cmap(normed_color(3)))
    ax.text(0.05, 0.14, r'7495-91xx', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes, color=cmap(normed_color(4)))
    ax.text(0.05, 0.08, r'7495-127xx', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes, color=cmap(normed_color(5)))

    ax.plot([0.8,0.8], ylim, dashes=[10,5], lw=2, color='k', zorder=2)

    lim = [-5, fgoodpix_saved.shape[0]+4]

    ax = fig.add_axes([0.45,0.13,0.22,0.22], facecolor='0.9')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.tick_params(which='both', bottom=False, left=False)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    cax = fig.add_axes([0.675, 0.13, 0.01, 0.22])
    img = ax.imshow(fgoodpix_saved.T, origin='lower', interpolation='nearest', cmap='gray',
                    vmin=0.5, vmax=None) #0.9)#,
                    #norm=colors.LogNorm(vmin=0.1, vmax=1.0))
    cb = pyplot.colorbar(img, cax=cax)
#    cb.locator = ticker.LogLocator(base=10, subs=(1.,2.,4.,))
#    cb.update_ticks()
    cax.tick_params(labelsize=10)
    cax.text(3.5, 0.95, r'$\delta\Lambda$', horizontalalignment='center',
             verticalalignment='center', transform=cax.transAxes, fontsize=12)

    ax.text(-0.06, 0.05, r'7495-12704', horizontalalignment='center', verticalalignment='bottom',
            transform=ax.transAxes, rotation='vertical')
    
    ofile='../ms/rev1/good_spaxel_growth_v2.pdf'
#    ofile=None
    if ofile is None:
        pyplot.show()
    else:
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



