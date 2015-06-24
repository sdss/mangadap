#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock

from mangadap.drpfile import drpfile
from mangadap.proc.TemplateLibrary import TemplateLibrary
from matplotlib import pyplot
from mangadap.util.par import TemplateLibraryParSet
from mangadap.util.defaults import default_dap_source

from mangadap.util.fileio import readfits_1dspec
from mangadap.util.instrument import spectral_resolution, spectrum_velocity_scale

#-----------------------------------------------------------------------------

def overplot_test():
    drpf = drpfile(7495, 12703, 'CUBE')

    tpl_lib = TemplateLibrary('M11-MILES', process=False)
    unmasked = numpy.invert( tpl_lib.tplbm.flagged(tpl_lib.hdu['MASK'].data) )

    print(tpl_lib.file_list)
    print(tpl_lib.ntpl)

    i = 0

    wave_0, flux_0 = readfits_1dspec(tpl_lib.file_list[i])

    wave_t = tpl_lib.hdu['WAVE'].data[i,unmasked[i,:]]
    flux_t = tpl_lib.hdu['FLUX'].data[i,unmasked[i,:]]

    print(wave_0.shape)
    print(wave_t.shape)
    print(flux_0.shape)
    print(flux_t.shape)

    pyplot.plot(tpl_lib.hdu['WAVE'].data[i,unmasked[i,:]],tpl_lib.hdu['FLUX'].data[i,unmasked[i,:]])
    pyplot.plot(wave_0, flux_0)
    pyplot.show()

    pyplot.plot(wave_0, flux_0-flux_t)
    pyplot.show()


def run_test():

#    drpf = drpfile(7495, 12703, 'CUBE')
    drpf = drpfile(7495, 12703, 'CUBE', read=True)

    sres = spectral_resolution(drpf.hdu['WAVE'].data, drpf.hdu['SPECRES'].data, log10=True)
    velscale = spectrum_velocity_scale(drpf.hdu['WAVE'].data, log10=True)

    print(drpf.file_path())

#    tpl_lib = TemplateLibrary('M11-MILES', drpf=drpf, directory_path='.', force=True)
    tpl_lib = TemplateLibrary('STELIB', velscale=velscale, sres=sres, directory_path='.',
                              processed_file='test.fits', force=True)
#    tpl_lib = TemplateLibrary('STELIB', process=False)

#    dapsrc = default_dap_source()
#
#    # M11-MARCS
#    template_libraries = TemplateLibraryParSet(key='M11-MARCS',
#                                file_search=dapsrc+'/external/templates/m11_marcs/*_s.fits',
#                                # TODO: This is the resolution in the header of the files, is it
#                                # right?
#                                fwhm=2.73, in_vacuum=False,
#                                wave_limit=numpy.array([ None, None ]),
#                                lower_flux_limit=None)
#    tpl_lib = TemplateLibrary('MILES-AVG', tpllib_list=template_libraries, process=False)
    
    print(tpl_lib.hdu['FLUX'].shape)
    nspec = tpl_lib.hdu['FLUX'].shape[0]

    unmasked = numpy.invert( tpl_lib.tplbm.flagged(tpl_lib.hdu['MASK'].data, flag=['NO_DATA','WAVE_INVALID','FLUX_INVALID']) )

#    edges = numpy.linspace(0,985,num=11,dtype=int)
#    for j in range(0,10):
#        for i in range(edges[j],edges[j+1]):
#        print(tpl_lib.file_list[i])

    for i in range(0,nspec):
#    for i in range(969,970):
#        print(tpl_lib.file_list[i])
#        pyplot.plot(tpl_lib.hdu['WAVE'].data[i,unmasked[i,:]], tpl_lib.hdu['FLUX'].data[i,unmasked[i,:]])
        x = tpl_lib.hdu['WAVE'].data[:]
#        x = tpl_lib.hdu['WAVE'].data[i,:]
        y = tpl_lib.hdu['FLUX'].data[i,:]
        pyplot.plot(x[unmasked[i,:]], y[unmasked[i,:]], ms=5, linestyle='none', marker='.', markeredgecolor=None, color='black')
        pyplot.plot(x[numpy.invert(unmasked[i,:])], y[numpy.invert(unmasked[i,:])], ms=10, linestyle='none', marker='.', markeredgecolor=None, color='red')
        df = (numpy.amax(y)-numpy.amin(y))*1.2
        f0 = (numpy.amax(y)+numpy.amin(y)-df)/2.0
        pyplot.ylim([f0, f0+df])
#        pyplot.plot(tpl_lib.hdu['WAVE'].data[unmasked[i,:]], tpl_lib.hdu['FLUX'].data[i,unmasked[i,:]])
        pyplot.show()

#    tpl_lib.hdu.info()
#    print(tpl_lib.hdu['WAVE'].shape)
#    print(tpl_lib.hdu['FLUX'].shape)
#    pyplot.plot(tpl_lib.hdu['WAVE'].data[101,:], tpl_lib.hdu['FLUX'].data[101,:])
#    pyplot.show()


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    #overplot_test()
    run_test()

    print('Elapsed time: {0} seconds'.format(clock() - t))



