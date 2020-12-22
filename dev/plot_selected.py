
import os
import time

from IPython import embed

import numpy

from astropy.io import fits

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot, rc, colors, ticker

from mangadap.util.drpfits import DRPFitsBitMask
from mangadap.dapfits import DAPCubeBitMask
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

def main():
    select_file = 'repr/representative_spectra_selection_v2.fits'
    spec_file = 'repr/benchmark_spectra_v2.fits'
    flag_file = 'repr/representative_spectra_flags_v2.db'

#    mod_file = None
#    plot_file = 'repr/representative_spectra_v2.pdf'

#    mod_file = 'repr/benchmark_spectra_v2_model.fits'
#    plot_file = 'repr/representative_spectra_v2_model.pdf'

#    mod_file = 'repr/run5/run5_MILESHC_4_8_MASTARSSP_1_14_ELBMPL9_ELPMPL11.fits.gz'
#    plot_file = 'repr/representative_spectra_v2_model_run5_ssp_14_sig1.4.pdf'
#    mod_file = 'repr/run5/run5_MILESHC_4_8_MASTARSSP_1_14_ELBMPL9_ELPMPL11A.fits.gz'
#    plot_file = 'repr/representative_spectra_v2_model_run5_ssp_14_sig3.0.pdf'
#    mod_file = 'repr/run5/run5_MILESHC_4_8_MASTARSSP_1_14_ELBMPL9_ELPMPL11B.fits.gz'
#    plot_file = 'repr/representative_spectra_v2_model_run5_ssp_14_sigmix.pdf'
    mod_file = 'repr/run7/run7_MILESHC_4_8_MASTARSSP_1_14_ELBMPL9_ELPMPL11UPA.fits.gz'
    plot_file = 'repr/representative_spectra_v2_model_run7_ssp_14_sigmixup.pdf'

    select_hdu = fits.open(select_file)
    spec_hdu = fits.open(spec_file)
    mod_hdu = None if mod_file is None else fits.open(mod_file)

    nspec = spec_hdu['FLUX'].shape[0]

    restwave_limits = numpy.array([3575, 10300])
    d4000_limits = numpy.array([3700, 4300])
    mg_limits = numpy.array([5100, 5285])
    halpha_limits = numpy.array([6460, 6676])
    ca_limits = numpy.array([8450, 8750])

    drp_bm = DRPFitsBitMask()
    cube_bm = DAPCubeBitMask()

    wave = spec_hdu['WAVE'].data
    flux = numpy.ma.MaskedArray(spec_hdu['FLUX'].data,
                                mask=drp_bm.flagged(spec_hdu['MASK'].data,
                                                flag=['DONOTUSE', 'FORESTAR']))

    if mod_hdu is None:
        model = numpy.ma.MaskedArray(spec_hdu['MODEL'].data,
                                     mask=cube_bm.flagged(spec_hdu['MODEL_MASK'].data,
                                                      flag='NOMODEL'))
        cont = numpy.ma.MaskedArray(spec_hdu['MODEL'].data - spec_hdu['EMLINE'].data,
                                    mask=cube_bm.flagged(spec_hdu['MODEL_MASK'].data,
                                                      flag='NOMODEL'))
        stellar = numpy.ma.MaskedArray(spec_hdu['STELLAR'].data,
                                       mask=cube_bm.flagged(spec_hdu['STELLAR_MASK'].data,
                                                      flag='NOMODEL'))
    else:
        model = numpy.ma.MaskedArray(mod_hdu['MODEL'].data, mask=mod_hdu['MODEL_MASK'].data > 0)
        model = StellarContinuumModel.reset_continuum_mask_window(model)
        cont = numpy.ma.MaskedArray(mod_hdu['MODEL'].data - mod_hdu['EMLINE'].data,
                                    mask=mod_hdu['MODEL_MASK'].data > 0)
        cont = StellarContinuumModel.reset_continuum_mask_window(cont)
        
        stellar = numpy.ma.MaskedArray(mod_hdu['STELLAR'].data,
                                       mask=mod_hdu['STELLAR_MASK'].data > 0)
        stellar = StellarContinuumModel.reset_continuum_mask_window(stellar)

    if not os.path.isfile(flag_file):
        include_flags = numpy.ones(nspec, dtype=int)
        numpy.savetxt(flag_file, numpy.array([select_hdu['PAR'].data['PLT'].astype(int),
                                              select_hdu['PAR'].data['IFU'].astype(int),
                                              select_hdu['PAR'].data['BIN'].astype(int),
                                              numpy.ones(nspec, dtype=int)]).T,
                      fmt=['%5d', '%5d', '%5d', '%3d'],
                      header='{0:>3} {1:>5} {2:>5} {3:>3}'.format('PLT', 'IFU', 'BIN', 'SEL'))

    rc('font', size=8)

    with PdfPages(plot_file) as pdf:

#        j = 0
        for i in range(nspec):
            print('{0}/{1}'.format((i+1),nspec), end='\r')
            lambda_limits = restwave_limits*(1+select_hdu['PAR'].data['Z'][i])
            indx = numpy.logical_not(numpy.ma.getmaskarray(flux)[i,:]) \
                        & (wave > lambda_limits[0]) & (wave < lambda_limits[1])
            if numpy.sum(indx) == 0:
                continue
            if numpy.sum(numpy.logical_not(numpy.ma.getmaskarray(model)[i,:])) == 0:
                continue

            mod_lim = [numpy.ma.amin(model[i,indx]), numpy.ma.amax(model[i,indx])]
            Df = (mod_lim[1]-mod_lim[0])*1.5
            flux_limits = numpy.mean(mod_lim) + numpy.array([-Df/2, Df/2])

            if not numpy.all(numpy.isfinite(flux_limits)):
                embed()
                exit()

            fig = pyplot.figure()

            # Full spectrum
            ax = fig.add_axes([0.06, 0.45, 0.92, 0.5])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.step(wave, flux[i,:], zorder=2, color='k', lw=0.5, where='mid')
            ax.plot(wave, stellar[i,:], zorder=3, color='C1', lw=1.0)
            ax.plot(wave, model[i,:], zorder=4, color='C3', lw=1.0)
            ax.plot(wave, cont[i,:], zorder=5, color='C9', lw=1.0)

            ax.text(-0.05, 0.5, r'$F_\lambda$',
                    ha='center', va='center', transform=ax.transAxes, rotation='vertical')
            ax.text(0.5, -0.1, r'$\lambda_{\rm obs}\ [{\rm \AA}]$',
                    ha='center', va='center', transform=ax.transAxes)

            ax.text(0.01, 1.03,
                    '{0}: {1}-{2}-{3}; '.format(i, select_hdu['PAR'].data['PLT'][i],
                                                select_hdu['PAR'].data['IFU'][i],
                                                select_hdu['PAR'].data['BIN'][i])
                    + 'z={0:.4f}; '.format(select_hdu['PAR'].data['Z'][i])
                    + 'S/N={0:.1f}; '.format(select_hdu['PAR'].data['SNR'][i])
                    + r'$\sigma_{\ast,{\rm obs}}$='
                    + r'{0:.1f}; D4000={1:.1f}; H$\alpha$ EW={2:.1f}'.format(
                        select_hdu['PAR'].data['SIGMA'][i], select_hdu['PAR'].data['D4000'][i],
                        select_hdu['PAR'].data['HAEW'][i]),
                    horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


            lambda_limits = d4000_limits*(1+select_hdu['PAR'].data['Z'][i])
            indx = numpy.logical_not(numpy.ma.getmaskarray(flux)[i,:]) \
                        & (wave > lambda_limits[0]) & (wave < lambda_limits[1])
            mod_lim = [0,1] if numpy.sum(indx) == 0 \
                        else [numpy.ma.amin(model[i,indx]), numpy.ma.amax(model[i,indx])]
            Df = (mod_lim[1]-mod_lim[0])*1.5
            flux_limits = numpy.mean(mod_lim) + numpy.array([-Df/2, Df/2])

            ax = fig.add_axes([0.06, 0.08, 0.19, 0.27])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.plot(wave, flux[i,:], zorder=2, color='k', lw=0.5)
            ax.plot(wave, stellar[i,:], zorder=3, color='C1', lw=1.0)
            ax.plot(wave, model[i,:], zorder=4, color='C3', lw=1.0)
            ax.plot(wave, cont[i,:], zorder=5, color='C9', lw=1.0)

            ax.text(0.5, -0.2, r'$\lambda_{\rm obs}\ [{\rm \AA}]$',
                    ha='center', va='center', transform=ax.transAxes)




            lambda_limits = mg_limits*(1+select_hdu['PAR'].data['Z'][i])
            indx = numpy.logical_not(numpy.ma.getmaskarray(flux)[i,:]) \
                        & (wave > lambda_limits[0]) & (wave < lambda_limits[1])
            mod_lim = [0,1] if numpy.sum(indx) == 0 \
                        else [numpy.ma.amin(model[i,indx]), numpy.ma.amax(model[i,indx])]
            Df = (mod_lim[1]-mod_lim[0])*1.5
            flux_limits = numpy.mean(mod_lim) + numpy.array([-Df/2, Df/2])

            ax = fig.add_axes([0.30, 0.08, 0.19, 0.27])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.plot(wave, flux[i,:], zorder=2, color='k', lw=0.5)
            ax.plot(wave, stellar[i,:], zorder=3, color='C1', lw=1.0)
            ax.plot(wave, model[i,:], zorder=4, color='C3', lw=1.0)
            ax.plot(wave, cont[i,:], zorder=5, color='C9', lw=1.0)

            ax.text(0.5, -0.2, r'$\lambda_{\rm obs}\ [{\rm \AA}]$',
                    ha='center', va='center', transform=ax.transAxes)




            lambda_limits = halpha_limits*(1+select_hdu['PAR'].data['Z'][i])
            indx = numpy.logical_not(numpy.ma.getmaskarray(flux)[i,:]) \
                        & (wave > lambda_limits[0]) & (wave < lambda_limits[1])
            mod_lim = [0,1] if numpy.sum(indx) == 0 \
                        else [numpy.ma.amin(model[i,indx]), numpy.ma.amax(model[i,indx])]
            Df = (mod_lim[1]-mod_lim[0])*1.5
            flux_limits = numpy.mean(mod_lim) + numpy.array([-Df/2, Df/2])

            ax = fig.add_axes([0.54, 0.08, 0.19, 0.27])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.plot(wave, flux[i,:], zorder=2, color='k', lw=0.5)
            ax.plot(wave, stellar[i,:], zorder=3, color='C1', lw=1.0)
            ax.plot(wave, model[i,:], zorder=4, color='C3', lw=1.0)
            ax.plot(wave, cont[i,:], zorder=5, color='C9', lw=1.0)

            ax.text(0.5, -0.2, r'$\lambda_{\rm obs}\ [{\rm \AA}]$',
                    ha='center', va='center', transform=ax.transAxes)




            lambda_limits = ca_limits*(1+select_hdu['PAR'].data['Z'][i])
            indx = numpy.logical_not(numpy.ma.getmaskarray(flux)[i,:]) \
                        & (wave > lambda_limits[0]) & (wave < lambda_limits[1])
            mod_lim = [0,1] if numpy.sum(indx) == 0 \
                        else [numpy.ma.amin(model[i,indx]), numpy.ma.amax(model[i,indx])]
            Df = (mod_lim[1]-mod_lim[0])*1.5
            flux_limits = numpy.mean(mod_lim) + numpy.array([-Df/2, Df/2])

            ax = fig.add_axes([0.78, 0.08, 0.19, 0.27])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(lambda_limits)
            ax.set_ylim(flux_limits)
            ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
            ax.plot(wave, flux[i,:], zorder=2, color='k', lw=0.5)
            ax.plot(wave, stellar[i,:], zorder=3, color='C1', lw=1.0)
            ax.plot(wave, model[i,:], zorder=4, color='C3', lw=1.0)
            ax.plot(wave, cont[i,:], zorder=5, color='C9', lw=1.0)

            ax.text(0.5, -0.2, r'$\lambda_{\rm obs}\ [{\rm \AA}]$',
                    ha='center', va='center', transform=ax.transAxes)


            pdf.savefig(orientation='landscape')
            fig.clear()
            pyplot.close(fig)
            
#            j += 1
#            if j == 10:
#                break

    print('{0}/{0}'.format(nspec))

if __name__ == '__main__':
    t = time.perf_counter()
    main()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



