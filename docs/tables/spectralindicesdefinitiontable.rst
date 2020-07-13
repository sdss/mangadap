=========  ===========  =========================================================================================
Key        Type         Description                                                                              
=========  ===========  =========================================================================================
TYPE       ``str_``     Type of spectral index, either absorption or bandhead.                                   
ID         ``int64``    ID number for the spectral index                                                         
NAME       ``str_``     Unique name for the spectral index                                                       
PASSBAND   ``float64``  Lower and upper rest wavelength of the main index passband (absorption-line indices only)
BLUEBAND   ``float64``  Lower and upper rest wavelength of the blue sideband used to define the linear continuum 
REDBAND    ``float64``  Lower and upper rest wavelength of the red sideband used to define the linear continuum  
UNIT       ``str_``     Index units                                                                              
COMPONENT  ``bool_``    Flag if the index is a component of a multi-component index (ignored)                    
INTEGRAND  ``str_``     The flux integrand for the index.                                                        
ORDER      ``str_``     The numerator-denominator order of the index (bandhead/color indices only)               
=========  ===========  =========================================================================================

