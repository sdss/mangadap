==============  ===========  ====================================================================================================================================================
Key             Type         Description                                                                                                                                         
==============  ===========  ====================================================================================================================================================
BINID           ``int64``    Bin ID number                                                                                                                                       
BINID_INDEX     ``int64``    0-indexed number of bin                                                                                                                             
MASK            ``int16``    Maskbit value                                                                                                                                       
BEGPIX          ``int64``    Index of the first pixel included in the fit                                                                                                        
ENDPIX          ``int64``    Index of the pixel just beyond the last pixel included in fit                                                                                       
NPIXTOT         ``int64``    Total number of pixels in the spectrum to be fit.                                                                                                   
NPIXFIT         ``int64``    Number of pixels used by the fit.                                                                                                                   
TPLWGT          ``float64``  Optimal weight of each template.                                                                                                                    
TPLWGTERR       ``float64``  Nominal error in the weight of each template.                                                                                                       
USETPL          ``bool_``    Flag that each template was included in the fit.                                                                                                    
ADDCOEF         ``float64``  Coefficients of the additive polynomial, if included.                                                                                               
MULTCOEF        ``float64``  Coefficients of the multiplicative polynomial, if included.                                                                                         
KININP          ``float64``  Input guesses for the kinematics                                                                                                                    
KIN             ``float64``  Best-fitting stellar kinematics                                                                                                                     
KINERR          ``float64``  Errors in the best-fitting stellar kinematics                                                                                                       
CHI2            ``float64``  Chi-square figure-of-merit for the fit                                                                                                              
RCHI2           ``float64``  Reduced chi-square figure-of-merit for the fit                                                                                                      
CHIGRW          ``float64``  Value of the error-normalized residuals at 0, 68%, 95%, 99%, and 100% growth                                                                        
RMS             ``float64``  Root-mean-square of the fit residuals.                                                                                                              
RMSGRW          ``float64``  Value of absolute value of the fit residuals at 0, 68%, 95%, 99%, and 100% growth                                                                   
FRMS            ``float64``  Root-mean-square of the fractional residuals (i.e., residuals/model).                                                                               
FRMSGRW         ``float64``  Value of absolute value of the fractional residuals at 0, 68%, 95%, 99%, 100% growth                                                                
SIGMACORR_SRES  ``float64``  Quadrature correction for the stellar velocity dispersion determined by the mean difference in spectral resolution between galaxy and template data.
SIGMACORR_EMP   ``float64``  Quadrature correciton for the stellar velocity dispersion determined by fitting the optimal template to one resolution matched to the galaxy data.  
==============  ===========  ====================================================================================================================================================

