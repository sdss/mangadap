
import numpy
from astropy.io import fits
from matplotlib import pyplot

from mangadap.util.bitmask import BitMask
from mangadap.config import defaults

sdssbits_file = defaults.sdss_maskbits_file()

targ1bm = BitMask.from_par_file(sdssbits_file, 'MANGA_TARGET1')
targ2bm = BitMask.from_par_file(sdssbits_file, 'MANGA_TARGET2')
targ3bm = BitMask.from_par_file(sdssbits_file, 'MANGA_TARGET3')

hdu = fits.open('drpall-v3_1_1.fits')

#print(hdu['MANGA'].columns.names)

indx = hdu['MANGA'].data['z'] < 0

mngtarg1 = hdu['MANGA'].data['mngtarg1'][indx]
mngtarg2 = hdu['MANGA'].data['mngtarg2'][indx]
mngtarg3 = hdu['MANGA'].data['mngtarg3'][indx]

print('MANGA_TARGET1')
for b in numpy.unique(mngtarg1):
    print(b, targ1bm.flagged_bits(b))

print('MANGA_TARGET2')
for b in numpy.unique(mngtarg2):
    print(b, targ2bm.flagged_bits(b))

print('MANGA_TARGET3')
for b in numpy.unique(mngtarg3):
    print(b, targ3bm.flagged_bits(b))

indx = numpy.where(indx)[0]
for i in indx:
    print('{0:>5} {1:>5} {2:12.5f} {3:8.2f} {4:7.1f} {5:7.1f} {6:>15} {7:>15} {8:>15}'.format(
          hdu['MANGA'].data['plate'][i], hdu['MANGA'].data['ifudsgn'][i], 
          hdu['MANGA'].data['z'][i], hdu['MANGA'].data['nsa_elpetro_ba'][i], 
          hdu['MANGA'].data['nsa_elpetro_phi'][i], hdu['MANGA'].data['nsa_elpetro_th50_r'][i], 
          ','.join(targ1bm.flagged_bits(hdu['MANGA'].data['mngtarg1'][i])),
          ','.join(targ2bm.flagged_bits(hdu['MANGA'].data['mngtarg2'][i])),
          ','.join(targ3bm.flagged_bits(hdu['MANGA'].data['mngtarg3'][i]))))

print(numpy.sum(indx))
exit()


indx = numpy.where(hdu['MANGA'].data['PLATEIFU'] == '8261-12705')[0][0]

print(hdu['MANGA'].data['z'][indx])
print(hdu['MANGA'].data['nsa_sersic_ba'][indx], hdu['MANGA'].data['nsa_elpetro_ba'][indx])
print(hdu['MANGA'].data['nsa_sersic_phi'][indx], hdu['MANGA'].data['nsa_elpetro_phi'][indx])
print(hdu['MANGA'].data['nsa_sersic_th50'][indx], hdu['MANGA'].data['nsa_elpetro_th50_r'][indx])
exit()

pyplot.scatter(hdu['MANGA'].data['nsa_sersic_ba'], hdu['MANGA'].data['nsa_elpetro_ba'],
               marker='.', s=30, lw=0, color='k')
pyplot.show()


