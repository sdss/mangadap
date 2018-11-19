# Galfit utils
from astropy.io import fits
import numpy as np
import numpy.ma as ma

def retrieve_PSF(plate, ifu, vor):
    #For now, just median stack for white-light image. Play aroubd with later
    hdulist3 = fits.open('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(ifu)+'-LOGCUBE.fits.gz')
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
    hdu.writeto('/Users/ppzaf/Work/MaNGA/dap/master/examples/%s-%s-input.fits' % (plate, ifu), clobber=True)
    #Save as input fits file
    ivars_masked = np.ma.array(ivars, mask=np.isclose(ivars, 0)) #Mask is so don't have to take sqrt(0)
    vars_masked = 1 / ivars_masked
    var_image = np.sum(vars_masked, axis=2)
    sig_image = np.sqrt(var_image)
    sig_image2 = np.ma.filled(sig_image, fill_value=1e9) #Can't save masked arrays as fits images...
    hdu1 = fits.PrimaryHDU(sig_image2, header=header) #Copy maps header into this coz need EXPTIME keyword.
    hdu1.writeto('/Users/ppzaf/Work/MaNGA/dap/master/examples/%s-%s-sigma.fits' % (plate, ifu), clobber=True)
    #Retrieve a PSF image  - for now, it's the median stacking of the four psf images, but this might change
    PSF_whitelight_im, grizpsf = retrieve_PSF(plate, ifu, vor)
    hdu2 = fits.PrimaryHDU(PSF_whitelight_im, header=header) #Copy maps header into this coz need EXPTIME keyword.
    hdu2.writeto('/Users/ppzaf/Work/MaNGA/dap/master/examples/%s-%s-PSF.fits' % (plate, ifu), clobber=True)
    #Make it write a feedme file
    mask_galfit = np.median(mask, axis=2) #badpixels>0, goodpixels==0. NOTE: Check this works always...
    hdu3 = fits.PrimaryHDU(mask_galfit, header=header) #Copy maps header into this coz need EXPTIME keyword.
    hdu3.writeto('/Users/ppzaf/Work/MaNGA/dap/master/examples/%s-%s-mask.fits' % (plate, ifu), clobber=True)
    return summed_image

def write_feedme(plate, ifu, summed_image, Rb, Rd, header):
    #Write a Galfit .feedme file for input into GalfitM
    Rbpix = 2*Rb
    Rdpix = 2*Rd
    f = open('/Users/ppzaf/Work/MaNGA/dap/master/examples/example-feedme', 'r')
    linelist = f.readlines()
    #Replace lines - atm for just an n=4 bulge and n=1 disk.
    linelist[2] = 'A) %s-%s-input.fits # Galaxy input image\n' % (plate, ifu)
    linelist[3] = 'B) %s-%s-output.fits #Output image block\n' % (plate, ifu)
    linelist[4] = 'C) %s-%s-sigma.fits # Galaxy sigma image\n' % (plate, ifu)
    linelist[5] = 'D) %s-%s-PSF.fits # PSF image\n' % (plate, ifu)
    linelist[7] = 'F) %s-%s-mask.fits	#Bad pixel mask\n' % (plate, ifu)
    linelist[9] = 'H) 0 %s 0 %s #Area to fit over\n' % (summed_image.shape[0], summed_image.shape[1])
    linelist[43] = '1) %s %s 1 1  #X and Y posn of gal (pix)\n' % (np.round(summed_image.shape[0]/2), np.round(summed_image.shape[1]/2))
    linelist[45] = '4) %s       1       #     R_e [Pixels]\n' % Rbpix #Was in arcsec, convert to pix
    linelist[47] = '9) %s       1       # axis ratio (b/a)\n' % (1-header['ECOOELL'])
    linelist[48] = '10) %s     1        # Position angle (deg)\n' % header['ECOOPA']
    linelist[54] = '1) %s %s 1 1 # X and Y posn of gal(pix)\n' % (np.round(summed_image.shape[0]/2), np.round(summed_image.shape[1]/2))
    linelist[56] = ' 4) %s       1       #     Rd [Pixels]\n' % Rdpix
    linelist[58] = '9) %s 1             #axis ratio of disk (b/a)\n' % (1-header['ECOOELL']) #Keep bulge and disk bas and PAs the same for now
    linelist[59] = '10) %s 1        #Position angle (deg)\n' % header['ECOOPA']
    with open('/Users/ppzaf/Work/MaNGA/dap/master/examples/%s-%s-feedme' % (plate, ifu), 'w') as file_handler:
        for item in linelist:
            file_handler.write("{}".format(item))
