from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import os
from os.path import join
import copy
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

from imp import reload

from mangadap import dap_access
from mangadap.plot import util
from mangadap.stack import select
from mangadap.stack import stack


plateifus = ['7443-1901', '7443-9101', '7443-12701']

"""SAMPLE SELECTION"""

out_bin = []
out = []
for plateifu in plateifus:
    # Specify File Name
    # filename = 'manga-7443-1901-LOGCUBE_BIN-RADIAL-004.fits'
    filename = 'manga-{}-LOGCUBE_BIN-RADIAL-004.fits'.format(plateifu)
    file_kws = util.parse_fits_filename(filename)
    path_data = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), os.getenv('MANGADRP_VER'),
                     os.getenv('MANGADAP_VER'), file_kws['plate'],
                     file_kws['ifudesign'])
    
    # Read in data
    gal = dap_access.DAPAccess(path_data, file_kws)
    gal.get_all_ext()
    
    """BIN SELECTION"""
    from mangadap.stack import select
    reload(select)
    # select bins with luminosity-weighted bin radius > 1 Re
    outer = gal.bins.binr > 1.0
    D4000_notnan = select.notnan(gal.sindx.indx.D4000, nanvals=-9999)
    D4000_sample = select.join_logical_and([outer, D4000_notnan])
    
    """STACKING"""
    from mangadap.stack import stack
    reload(stack)
    outer_Ha = gal.flux_ew.Ha6564.loc[outer].mean()
    outer_D4000 = gal.sindx.indx.D4000.loc[D4000_sample].mean()
    
    out.append([gal.mangaid, gal.header['OBJRA'], gal.header['OBJDEC'],
               outer_Ha, outer_D4000])

    out_bin.append([gal.flux_ew.Ha6564.loc[outer],
                   gal.sindx.indx.D4000.loc[D4000_sample]])

"""Galaxy-internal stacking"""
# Do this in a function
df = pd.DataFrame(out, columns=['mangaid', 'RA', 'DEC', 'Ha', 'D4000'],
                  index=plateifus)

"""Cross-sample stacking"""
# Do the next three lines in a function
out_bin_T = [list(it) for it in zip(*out_bin)] # "transpose" list
out_bin_concat = [pd.concat(it) for it in out_bin_T] # concat data by type
df_bin = pd.concat(out_bin_concat, axis=1, keys=['Ha', 'D4000'])

# Print results
df_bin.Ha.mean()
df_bin.D4000.mean()





"""
Stacking options:
For now, assume dataframes and use built-in mean() and median() if possible.
"""




"""Another example"""
# user-defined set of bins
ind_bins = np.array([0, 1, 100, 200, 300, 350])
# convert to boolean index array
bins = select.int_to_bool_index(ind_bins, gal.flux_ew.Ha6564.shape)
# bin Halpha flux selection cut
high_halpha = gal.flux_ew.Ha6564 > gal.flux_ew.Ha6564.median()
# list of conditions (including remove bins where Halpha flux is NaN)
conditions = [bins, high_halpha, select.notnan(gal.flux_ew.Ha6564)]
# join conditions
sample = select.join_logical_and(conditions)


"""Ivar wtmean Ha flux and D4000 for > 1 Re"""
outer_Ha = stack.ivar_wtmean(gal.flux_ew.Ha6564, gal.fluxerr_ew.Ha6564, outer)

