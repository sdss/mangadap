=============  ===  =====================================================================================================================================================
Key            Bit  Description                                                                                                                                          
=============  ===  =====================================================================================================================================================
DIDNOTUSE      0    Pixel was ignored because it was flagged as either DONOTUSE or FORESTAR by the data reduction pipeline (DRP).                                        
FORESTAR       1    Pixel was ignored because it was flagged as FORESTAR by the data reduction pipeline (DRP).                                                           
LOW_SPECCOV    2    Pixel was ignored because the fraction of valid spectral channels limited the spectral coverage below the set threshold; see header keyword FSPECCOV.
LOW_SNR        3    Pixel was ignored because the S/N estimate of the spectrum was below the set threshold; see header keyword BINMINSN.                                 
NONE_IN_STACK  4    No valid pixels were available or valid for the stacked spectrum.                                                                                    
NO_STDDEV      5    Insufficient pixels to calculate the standard deviation in the stacked spectrum.                                                                     
IVARINVALID    6    Pixel ignored because inverse variance invalid.                                                                                                      
=============  ===  =====================================================================================================================================================

