#!/usr/bin/env python3

import os
import time
import numpy

from astropy.io import fits
from matplotlib import pyplot

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()
    print('Elapsed time: {0} seconds'.format(time.clock() - t))

    calls = [ './fit_selected_spectra.py benchmark_spectra_d4000_sig_haew_snr.fits.gz fit_dshn_mn-1.fits.gz --ppxf_iteration no_global_wrej --spec_flags d4000_sigma_haew_snr_include.db --degree -1 ' ]

    orders = numpy.arange(29)
    fboxes = 3127//(orders+1)

    for o,f in zip(orders, fboxes):

        calls += [ './fit_selected_spectra.py benchmark_spectra_d4000_sig_haew_snr.fits.gz fit_dshn_mn{0}.fits.gz --ppxf_iteration no_global_wrej --spec_flags d4000_sigma_haew_snr_include.db --degree {1}'.format(str(o).zfill(2),o) ]

        calls += [ './fit_selected_spectra.py benchmark_spectra_d4000_sig_haew_snr.fits.gz fit_dshn_fs{0}.fits.gz --ppxf_iteration fit_reject_filter --spec_flags d4000_sigma_haew_snr_include.db --filt_degree -1 --filter 1 --filter_opt subtract --filter_boxcar {1}'.format(str(f).zfill(4),f) ]

    nthreads = 6
    nfiles = min(nthreads, len(calls))
    
    of = [ open('fit_selected_{0}.scr'.format(str(i+1).zfill(2)), 'w') for i in range(nfiles) ]
    for i in range(nfiles):
        of[i].write('\n')
   
    i = 0
    for c in calls:
        of[i].write(c)
        of[i].write('\n')
        i = i + 1 if i < nfiles-1 else 0
    
    for i in range(nfiles):
        of[i].write('\n')
        of[i].close()
    


