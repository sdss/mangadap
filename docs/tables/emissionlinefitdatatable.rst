===========  ===========  ==========================================================================================================================================================
Key          Type         Description                                                                                                                                               
===========  ===========  ==========================================================================================================================================================
BINID        ``int64``    Bin ID number                                                                                                                                             
BINID_INDEX  ``int64``    0-indexed number of bin                                                                                                                                   
FIT_INDEX    ``int64``    The index in the fit database associated with each emission line.                                                                                         
MASK         ``int16``    Maskbit value for each emission line.                                                                                                                     
FLUX         ``float64``  The best-fitting flux of the emission line.                                                                                                               
FLUXERR      ``float64``  The error in the best-fitting emission-line flux                                                                                                          
KIN          ``float64``  The best-fitting kinematics in each emission line                                                                                                         
KINERR       ``float64``  The error in the best-fitting emission-line kinematics                                                                                                    
SIGMACORR    ``float64``  Quadrature correction in the emission-line velocity dispersion                                                                                            
SIGMAINST    ``float64``  Dispersion of the instrumental line-spread function at the location of each emission line.                                                                
SIGMATPL     ``float64``  Dispersion of the instrumental line-spread function of the emission-line templates.                                                                       
CONTAPLY     ``float64``  The value of any additive polynomial included in the fit at the location of each emission line                                                            
CONTMPLY     ``float64``  The value of any multiplicative polynomial included in the fit at the location of each emission line                                                      
CONTRFIT     ``float64``  The value of any extinction curve included in the fit at the location of each emission line                                                               
LINE_PIXC    ``int64``    The integer pixel nearest the center of each emission line.                                                                                               
AMP          ``float64``  The best-fitting amplitude of the emission line.                                                                                                          
ANR          ``float64``  The amplitude-to-noise ratio defined as the model amplitude divided by the median noise in the two (blue and red) sidebands defined for the emission line.
LINE_NSTAT   ``int64``    The number of pixels included in the fit metric calculations (LINE_RMS, LINE_FRMS, LINE_CHI2) near each emission line.                                    
LINE_RMS     ``float64``  The root-mean-square residual of the model fit near each emission line.                                                                                   
LINE_FRMS    ``float64``  The root-mean-square of the fractional residuals of the model fit near each emission line.                                                                
LINE_CHI2    ``float64``  The chi-square of the model fit near each emission line.                                                                                                  
BMED         ``float64``  The median flux in the blue sideband of each emission line                                                                                                
RMED         ``float64``  The median flux in the red sideband of each emission line                                                                                                 
EWCONT       ``float64``  The continuum value interpolated at the emission-line center (in the observed frame) used for the equivalent width measurement.                           
EW           ``float64``  The equivalent width of each emission line                                                                                                                
EWERR        ``float64``  The error in the equivalent width of each emission line                                                                                                   
===========  ===========  ==========================================================================================================================================================

