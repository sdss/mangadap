from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

from imp import reload

from mangadap import dap_access
from mangadap.plot import util
from mangadap.stack import select

# Specify File Name
filename = 'manga-7443-1901-LOGCUBE_BIN-NONE-003.fits'
file_kws = util.parse_fits_filename(filename)
path_data = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), os.getenv('MANGADRP_VER'),
                 os.getenv('MANGADAP_VER'), file_kws['plate'],
                 file_kws['ifudesign'])

# Read in data
gal = dap_access.DAPAccess(path_data, file_kws)
gal.get_all_ext()




from mangadap.stack import select
reload(select)
bins = np.array([0, 1])
indbool_bins = select.int_to_bool_index(bins, gal.flux_ew.Ha6564.shape)


high_halpha = gal.flux_ew.Ha6564 > gal.flux_ew.Ha6564.median()
mask = select.join_conditions([high_halpha, select.notnan(gal.flux_ew.Ha6564)])
gal.flux_ew.OII3727.loc[high_halpha].mean()



# Stack DAP Values from One Galaxy
def stack(val, bins, method='mean'):
    """Stack DAP values."""
    out = np.__dict__[method](val[bins])
    return out
    

a = stack(val=gal.flux_ew.OII3727, bins=np.array([0, 1]))
print(a)
print(np.mean(gal.flux_ew.OII3727[:2]))


high_halpha = gal.flux_ew.Ha6564 > gal.flux_ew.Ha6564.median()
low_halpha = gal.flux_ew.Ha6564 <= gal.flux_ew.Ha6564.median()
print('High Halpha: F(OII3727) =', gal.flux_ew.OII3727.loc[high_halpha].mean())



# The problem with stacking spectral indices is that -9999 was used as a NaN.
# The specind errors are also -9999 for these bins.
ind = np.where((gal.sindx.indx == -9999))
#ind
#gal.sindx.indx.FeH0p99
gal.sindx.indxerr.values[ind]


# Masking

# Drop NaNs
notnan = (gal.sindx.indx.D4000 != -9999.)
ind_D4000 = np.logical_and(low_halpha, notnan)
gal.sindx.indx.D4000.loc[ind_D4000].mean()


# Weighting

# inverse variance weighting
gal.sindx.indx.D4000.loc[ind_D4000].dot(gal.sindx.indxerr.D4000[ind_D4000]**-2.) / np.sum(gal.sindx.indxerr.D4000[ind_D4000]**-2.)


# ## Plotting

# Quick Plot
# Make a simple plot from the columns from a DataFrame with [pandas plotting functions](http://pandas.pydata.org/pandas-docs/stable/visualization.html).
gal.sindx.indx.D4000[ind_D4000].plot()


D4000 = pd.DataFrame({'val': gal.sindx.indx.D4000[ind_D4000],
                      'err': gal.sindx.indxerr.D4000[ind_D4000]*10.})

D4000.plot(kind='scatter', x='val', y='val', yerr='err')


