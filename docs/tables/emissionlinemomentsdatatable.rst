===========  ===========  ===================================================================
Key          Type         Description                                                        
===========  ===========  ===================================================================
BINID        ``int64``    Spectrum/Bin ID number                                             
BINID_INDEX  ``int64``    Index of the spectrum in the list of provided spectra.             
REDSHIFT     ``float64``  Redshift used for shifting the passbands                           
MASK         ``bool_``    Bad-value boolean or bit mask value for the moments                
BCEN         ``float64``  Center of the blue sideband.                                       
BCONT        ``float64``  Pseudo-continuum in the blue sideband                              
BCONTERR     ``float64``  Error in the blue-sideband pseudo-continuum                        
RCEN         ``float64``  Center of the red sideband.                                        
RCONT        ``float64``  Pseudo-continuum in the red sideband                               
RCONTERR     ``float64``  Error in the red-sideband pseudo-continuum                         
CNTSLOPE     ``float64``  Continuum slope used to determine the continuum at the line center.
FLUX         ``float64``  Summed flux (0th moment)                                           
FLUXERR      ``float64``  Error in the summed flux                                           
MOM1         ``float64``  Line centroid redshift (:math:`cz`; 1st moment)                    
MOM1ERR      ``float64``  Error in the line centroid redshift                                
MOM2         ``float64``  Line standard deviation (2nd moment)                               
MOM2ERR      ``float64``  Error in the line standard deviation                               
SINST        ``float64``  Instrumental dispersion at the line centroid                       
BMED         ``float64``  Median flux in the blue sideband used for EW                       
RMED         ``float64``  Median flux in the red sideband used for EW                        
EWCONT       ``float64``  Continuum value used for EW calculation                            
EW           ``float64``  Equivalent width (FLUX/pseudo continuum)                           
EWERR        ``float64``  Error in the equivalent width                                      
===========  ===========  ===================================================================

