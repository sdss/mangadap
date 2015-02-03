file='M11_STELIB_res3.1_imf_ss_age_0.200000_z_0.01.fits'
data = readfits(file, hdr)
crpix = fxpar(hdr, 'CRPIX1')
crval = fxpar(hdr, 'CRVAL1')
cdelt = fxpar(hdr, 'CDELT1')
npix = n_elements(data)
wave = crval+(dindgen(npix)+1-crpix)*cdelt

file_s='M11_STELIB_res3.1_imf_ss_age_0.200000_z_0.01_s.fits'
data_s = readfits(file_s, hdr)
crpix = fxpar(hdr, 'CRPIX1')
crval = fxpar(hdr, 'CRVAL1')
cdelt = fxpar(hdr, 'CDELT1')
npix = n_elements(data_s)
wave_s = crval+(dindgen(npix)+1-crpix)*cdelt




