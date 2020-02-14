========================  ===  ======================================================================================================================================================================
Key                       Bit  Description                                                                                                                                                           
========================  ===  ======================================================================================================================================================================
DIDNOTUSE                 0    Pixel was ignored because it was flagged as DONOTUSE or FORESTAR by the DRP, or as LOW_SPECCOV, LOW_SNR, or NONE_IN_STACK in the binning step.                        
FORESTAR                  1    Pixel was ignored because it was flagged as FORESTAR by the DRP.                                                                                                      
LOW_SNR                   2    Pixel was ignored because the S/N estimate of the spectrum was below the set threshold; see header keyword ELMMINSN.                                                  
MAIN_EMPTY                5    Line moments not measured because the main passband had no observed fluxes.                                                                                           
BLUE_EMPTY                6    Line moments not measured because the blue passband had no observed fluxes.                                                                                           
RED_EMPTY                 7    Line moments not measured because the red passband had no observed fluxes.                                                                                            
MAIN_INCOMP               8    There were masked fluxes in the main passband.                                                                                                                        
BLUE_INCOMP               9    There were masked fluxes in the blue passband.                                                                                                                        
RED_INCOMP                10   There were masked fluxes in the red passband.                                                                                                                         
DIVBYZERO                 11   Could not compute emission-line moment or its error because of a zero devision.                                                                                       
NO_ABSORPTION_CORRECTION  12   Moment has not been corrected for underlying absorption because no stellar-continuum model was provided or the model did not cover the necessary spectral range.      
MAIN_JUMP                 13   Some fraction of, but less than all, the pixels in the main passband were not fit by the stellar-continuum model, meaning the bandpass integral should not be trusted.
BLUE_JUMP                 14   Some fraction of, but less than all, the pixels in the blue sideband were not fit by the stellar-continuum model, meaning the bandpass integral should not be trusted.
RED_JUMP                  15   Some fraction of, but less than all, the pixels in the red sideband were not fit by the stellar-continuum model, meaning the bandpass integral should not be trusted. 
JUMP_BTWN_SIDEBANDS       16   The stellar-continuum was subtracted from one sideband, but not the other, meaning that the determination of the continuum beneath the line is invalid.               
UNDEFINED_MOM1            17   First moment is undefined, likely because the passband was empty                                                                                                      
UNDEFINED_MOM2            18   Second moment is undefined, likely because it leads to the square-root of a negative number.                                                                          
UNDEFINED_BANDS           19   No moments calculated because the bandpasses were not, or improperly, defined.                                                                                        
NON_POSITIVE_CONTINUUM    20   Equivalent width measurements were not computed because the continuum in either the blue or red sidebands was not positive.                                           
========================  ===  ======================================================================================================================================================================

