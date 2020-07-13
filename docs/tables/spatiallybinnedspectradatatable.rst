==========  ===========  ===========================================================================================
Key         Type         Description                                                                                
==========  ===========  ===========================================================================================
BINID       ``int64``    Bin ID number                                                                              
NBIN        ``int64``    Number of spaxels in the bin                                                               
SKY_COO     ``float64``  The mean on-sky coordinates of the binned spaxels.                                         
LW_SKY_COO  ``float64``  The luminosity-weighted mean on-sky coordinates of the binned spaxels.                     
ELL_COO     ``float64``  The mean elliptical coordinates of the binned spaxels.                                     
LW_ELL_COO  ``float64``  The luminosity-weighted mean elliptical coordinates of the binned spaxels.                 
AREA        ``float64``  Total on-sky area of the bin                                                               
AREA_FRAC   ``float64``  Fraction of the expected area covered by good spaxels (not relevant to all binning schemes)
SIGNAL      ``float64``  Mean flux per pixel in the binned spectrum.                                                
VARIANCE    ``float64``  Mean variance per pixel in the binned spectrum.                                            
SNR         ``float64``  Mean S/N per pixel in the binned spectrum.                                                 
==========  ===========  ===========================================================================================

