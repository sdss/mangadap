# Galfit utils
from astropy.io import fits
import numpy as np
import numpy.ma as ma

def retrieve_PSF(plate, ifu, vor):
    #For now, just median stack for white-light image. Play aroubd with later
    hdulist3 = fits.open('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(ifu)+'-LOGCUBE.fits.gz')
    gpsf = hdulist3['GPSF'].data #was 16
    rpsf = hdulist3['RPSF'].data #was 17
    ipsf = hdulist3['IPSF'].data #was 18
    zpsf = hdulist3['ZPSF'].data #was 19
    grizpsf = np.array([gpsf,rpsf,ipsf,zpsf])
    psf_im = np.median(grizpsf, axis=0)
    #os.remove('./manga-'+str(plate)+'-'+str(ifu)+'-LOGCUBE-'+str(vor)+'-GAU-MILESHC.fits.gz')
    return psf_im, grizpsf

def broadband_images(fluxes, ivars, mask, header, plate, ifu, vor):
    #Collapse datacube into flat image for input into Galfit.
    summed_image = fluxes.sum(axis=2) #Must have an exptime image header
    hdu = fits.PrimaryHDU(summed_image, header=header) #Copy maps header into this coz beed EXPTIME keyword.
    hdu.writeto('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/%s-%s-input.fits' % (plate, ifu), clobber=True)
    #Save as input fits file
    ivars_masked = np.ma.array(ivars, mask=np.isclose(ivars, 0)) #Mask is so don't have to take sqrt(0)
    vars_masked = 1 / ivars_masked
    var_image = np.sum(vars_masked, axis=2)
    sig_image = np.sqrt(var_image)
    sig_image2 = np.ma.filled(sig_image, fill_value=1e9) #Can't save masked arrays as fits images...
    hdu1 = fits.PrimaryHDU(sig_image2, header=header) #Copy maps header into this coz need EXPTIME keyword.
    hdu1.writeto('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/%s-%s-sigma.fits' % (plate, ifu), clobber=True)
    #Retrieve a PSF image  - for now, it's the median stacking of the four psf images, but this might change
    PSF_whitelight_im, grizpsf = retrieve_PSF(plate, ifu, vor)
    hdu2 = fits.PrimaryHDU(PSF_whitelight_im, header=header) #Copy maps header into this coz need EXPTIME keyword.
    hdu2.writeto('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/%s-%s-PSF.fits' % (plate, ifu), clobber=True)
    #Make it write a feedme file
    mask_galfit = np.median(mask, axis=2) #badpixels>0, goodpixels==0. NOTE: Check this works always...
    hdu3 = fits.PrimaryHDU(mask_galfit, header=header) #Copy maps header into this coz need EXPTIME keyword.
    hdu3.writeto('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/%s-%s-mask.fits' % (plate, ifu), clobber=True)
