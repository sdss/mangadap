
import os
import numpy as np
import astropy.constants
from astropy.io import fits
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import matplotlib.colors as mcol

from mangadap.dapfits import DAPCubeBitMask
from mangadap.util.fileio import channel_dictionary
import marvin.utils.plot.map as mapplot
import pdb


plate = 7443
ifu = 1902
nsa_redshift = 0.0189259
vel = nsa_redshift * astropy.constants.c.to('km/s').value
plate_str = "{:.0f}".format(plate)
ifu_str = "{:.0f}".format(ifu)
outfil = 'check_new_binning_scheme_'+plate_str+'-'+ifu_str+'.pdf'

analysis_path = os.getcwd()+'/5x5n-GAU-MILESHC/'
fits_cube_fil = analysis_path+plate_str+'/'+ifu_str+'/manga-'+plate_str+'-'+ifu_str+'-LOGCUBE-5x5n-GAU-MILESHC.fits.gz'
fits_maps_fil = analysis_path+plate_str+'/'+ifu_str+'/manga-'+plate_str+'-'+ifu_str+'-MAPS-5x5n-GAU-MILESHC.fits.gz'
#pdb.set_trace()

# Read in CUBE
hdu_cube = fits.open(fits_cube_fil)
bm = DAPCubeBitMask()

wave = hdu_cube['WAVE'].data
flux = np.ma.MaskedArray(hdu_cube['FLUX'].data,
                         mask=bm.flagged(hdu_cube['MASK'].data,
                        ['IGNORED', 'FLUXINVALID', 'IVARINVALID', 'ARTIFACT']))
ivar = np.ma.MaskedArray(hdu_cube['IVAR'].data,
                         mask=bm.flagged(hdu_cube['MASK'].data,
                        ['IGNORED', 'FLUXINVALID', 'IVARINVALID', 'ARTIFACT']))
model = np.ma.MaskedArray(hdu_cube['MODEL'].data,
                          mask=bm.flagged(hdu_cube['MASK'].data, 'FITIGNORED'))

# Read in MAP
hdu_maps = fits.open(fits_maps_fil)
binid = hdu_maps['BINID'].data[1,:,:]
mask_ext = hdu_maps['EMLINE_GVEL'].header['QUALDATA']
emlc = channel_dictionary(hdu_maps, 'EMLINE_GVEL')
emission_vfield = np.ma.MaskedArray(hdu_maps['EMLINE_GVEL'].data[emlc['Ha-6564'],:,:],
                                    mask=hdu_maps[mask_ext].data[emlc['Ha-6564'],:,:] > 0)
emission_vfield_ivar = np.ma.MaskedArray(hdu_maps['EMLINE_GVEL_IVAR'].data[emlc['Ha-6564'],:,:],
                                        mask=hdu_maps[mask_ext].data[emlc['Ha-6564'],:,:] > 0)


mask_ext = hdu_maps['STELLAR_VEL'].header['QUALDATA']
stellar_vfield = np.ma.MaskedArray(hdu_maps['STELLAR_VEL'].data,
                                    mask=hdu_maps[mask_ext].data > 0)
stellar_vfield_ivar = np.ma.MaskedArray(hdu_maps['STELLAR_VEL_IVAR'].data,
                                        mask=hdu_maps[mask_ext].data > 0)


#Plot map of binids, emission line velocities, etc
fig, axes = pl.subplots(1, 3, figsize=(12,4))

rand_cmap = mcol.ListedColormap(np.random.rand(400,3))
mapplot.plot(value=binid, title='Bin ID', cmap=rand_cmap, fig=fig, ax=axes[0], cbrange=[0,315])

mapplot.plot(value=emission_vfield, title='Emission Line Velocity', cmap=cm.coolwarm,
             cbrange=[-120,120], fig=fig, ax=axes[1])

mapplot.plot(value=stellar_vfield, title='Stellar Velocity', cmap=cm.coolwarm,
             cbrange=[-120,120], fig=fig, ax=axes[2])

fig.tight_layout()
pl.savefig(outfil)
pdb.set_trace()