==============  ===  ====================================================================================================================================================================
Key             Bit  Description                                                                                                                                                         
==============  ===  ====================================================================================================================================================================
NO_DATA         0    Pixel has no data                                                                                                                                                   
WAVE_INVALID    1    Used to designate pixels in the 1D spectra that are outside the valid wavelength range defined by default_template_libraries().                                     
FLUX_INVALID    2    Used to designate pixels in the 1D spectra that are below the valid flux limit defined by default_template_libraries().                                             
SPECRES_EXTRAP  3    The spectral resolution has been matched to a value that was an extrapolation of the target spectral resolution samples.                                            
SPECRES_LOW     4    The spectral resolution was *not* matched to the target value because the target value was *higher* than the existing spectral resolution.                          
SPECRES_2PIXEL  5    Resampling has lead to the resolution being below the two pixel limit.  The spectral resolution has been changed to be exactly the two pixel limit for these pixels.
SPECRES_NOFLUX  6    Resolution matching has lead a template flux <= 0.                                                                                                                  
==============  ===  ====================================================================================================================================================================

