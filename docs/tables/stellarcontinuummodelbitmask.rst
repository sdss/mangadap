==================  ===  ===================================================================================================================================================================
Key                 Bit  Description                                                                                                                                                        
==================  ===  ===================================================================================================================================================================
DIDNOTUSE           0    Pixel was ignored because it was flagged in the binning step; see the DRPFits.do_not_stack_flags() and SpatiallyBinnedSpectra.do_not_fit_flags().                  
FORESTAR            1    Pixel was ignored because it was flagged as FORESTAR by the DRP.                                                                                                   
LOW_SNR             2    Pixel was ignored because the S/N estimate of the spectrum was below the set threshold; see header keyword SCMINSN.                                                
ARTIFACT            3    Pixel was ignored during the stellar-continuum fit because it was designated as containing an artifact.                                                            
OUTSIDE_RANGE       4    Pixel was ignored during the stellar-continuum fit because it was outside the designated spectral range of the fit; see the header keyword FITWAVE.                
EML_REGION          5    Pixel was ignored during the stellar-continuum fit because it contains an emission-line.                                                                           
TPL_PIXELS          6    These pixels were removed to ensure that the number of template spectral pixels was >= the number of fitted object pixels during the stellar-continuum fit.        
TRUNCATED           7    This region was truncated to avoid convolution errors at the edges of the spectral range during the stellar-continuum fit.                                         
PPXF_REJECT         8    Pixel was rejected during the pPXF fit via the clean parameter.                                                                                                    
INVALID_ERROR       9    Pixel ignored because the flux error was invalid.                                                                                                                  
NO_FIT              10   Spectrum not fit such that not parameters are provided.                                                                                                            
MAXITER             11   The fit optimizer reached the maximum number of iterations during the fit, which may imply failure to converge.                                                    
LARGE_CHI2          12   This pixel has a large chi^2 (definition of large TBD).                                                                                                            
LARGE_RESID         13   This pixel has a large residual (definition of large TBD).                                                                                                         
INSUFFICIENT_DATA   14   There were insufficient data in the fitting window to fit the line profile(s).                                                                                     
FIT_FAILED          15   Stellar-continuum fit failed according to status returned by scipy.optimize.least_squares.                                                                         
NEAR_BOUND          16   Fit parameters are within at or near the imposed boundaries in parameter space.                                                                                    
NEGATIVE_WEIGHTS    17   Weights for templates were negative.                                                                                                                               
BAD_SIGMA           18   Corrected velocity dispersion is below 50 km/s or above 400 km/s                                                                                                   
MIN_SIGMA           19   The fitted velocity dispersion is at the minimum allowed by pPXF (1/100th of a pixel)                                                                              
BAD_SIGMACORR_SRES  20   The instrumental dispersion correction for the velocity dispersion, based on the quadrature difference in the spectral resolution, is invalid.                     
BAD_SIGMACORR_EMP   21   The instrumental dispersion correction for the velocity dispersion, based on fitting the native resolution template to the resolution matched template, is invalid.
==================  ===  ===================================================================================================================================================================

