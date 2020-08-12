
import os
import time

import numpy

from astropy.io import fits

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot

from mangadap.util.drpfits import DRPFitsBitMask
from mangadap.dapfits import DAPCubeBitMask

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

def main():
    select_file = 'representative_spectra_selection.fits'
    spec_file = 'benchmark_spectra.fits'
    flag_file = 'representative_spectra_flags.db'
    plot_file = 'representative_spectra.pdf'

    select_hdu = fits.open(select_file)
    spec_hdu = fits.open(spec_file)

    nspec = spec_hdu['FLUX'].shape[0]

    restwave_limits = numpy.array([3575, 10300])

    drp_bm = DRPFitsBitMask()
    cube_bm = DAPCubeBitMask()

    wave = spec_hdu['WAVE'].data
    flux = numpy.ma.MaskedArray(spec_hdu['FLUX'].data,
                                mask=drp_bm.flagged(spec_hdu['MASK'].data,
                                                flag=['DONOTUSE', 'FORESTAR']))
    model = numpy.ma.MaskedArray(spec_hdu['MODEL'].data,
                                 mask=cube_bm.flagged(spec_hdu['MODEL_MASK'].data,
                                                      flag='NOMODEL'))
    cont = numpy.ma.MaskedArray(spec_hdu['MODEL'].data - spec_hdu['EMLINE'].data,
                                 mask=cube_bm.flagged(spec_hdu['MODEL_MASK'].data,
                                                      flag='NOMODEL'))

    if not os.path.isfile(flag_file):
        include_flags = numpy.ones(nspec, dtype=int)
        numpy.savetxt(flag_file, numpy.array([select_hdu['PAR'].data['PLT'].astype(int),
                                              select_hdu['PAR'].data['IFU'].astype(int),
                                              select_hdu['PAR'].data['BIN'].astype(int),
                                              numpy.ones(nspec, dtype=int)]).T,
                      fmt=['%5d', '%5d', '%5d', '%3d'],
                      header='{0:>3} {1:>5} {2:>5} {3:>3}'.format('PLT', 'IFU', 'BIN', 'SEL'))

    with PdfPages(plot_file) as pdf:

#        for i in range(nspec):
        for i in range(5):
            print('{0}/{1}'.format((i+1),nspec), end='\r')
            lambda_limits = restwave_limits*(1+select_hdu['PAR'].data['Z'][i])
            indx = numpy.invert(flux.mask[i,:]) & (spec_hdu['WAVE'].data > lambda_limits[0]) \
                        & (spec_hdu['WAVE'].data < lambda_limits[1])
            if numpy.sum(indx) == 0:
                continue
            srt=numpy.argsort(flux[i,indx])
            Df = (flux[i,indx][srt[int(0.9*numpy.sum(indx))]]
                    - flux[i,indx][srt[int(0.1*numpy.sum(indx))]])*3
            flux_limits = [ flux[i,indx][srt[int(0.5*numpy.sum(indx))]] - Df/2, 0]
            flux_limits[1] = flux_limits[0] + Df

            fig = pyplot.figure()

            ax = fig.add_axes([0.15, 0.3, 0.8, 0.6])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10)
            ax.tick_params(which='minor', length=5)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.plot(wave, flux[i,:], zorder=2, color='k', lw=0.5)
            ax.plot(wave, model[i,:], zorder=3, color='C3', lw=1.0)
            ax.plot(wave, cont[i,:], zorder=3, color='C0', lw=1.0)

            ax.set_xlabel(r'$\lambda$')
            ax.set_ylabel(r'$F(\lambda)$')

            ax.text(0.01, 1.03,
                    r'{0}: {1}-{2}-{3}; S/N={4:.1f}; $\sigma$={5:.1f}; D4000={6:.1f}; H$\alpha$ EW={7:.1f}'.format(
                            i+1, select_hdu['PAR'].data['PLT'][i], select_hdu['PAR'].data['IFU'][i],
                            select_hdu['PAR'].data['BIN'][i], select_hdu['PAR'].data['SNR'][i],
                            select_hdu['PAR'].data['SIGMA'][i], select_hdu['PAR'].data['D4000'][i],
                            select_hdu['PAR'].data['HAEW'][i]),
                    horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

            pdf.savefig(orientation='landscape')
            fig.clear()
            pyplot.close(fig)
        print('{0}/{0}'.format(nspec))

if __name__ == '__main__':
    t = time.perf_counter()
    main()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



