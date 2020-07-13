=========  ===========  ===============================================================================
Key        Type         Description                                                                    
=========  ===========  ===============================================================================
ID         ``int64``    Emission line ID number                                                        
NAME       ``str_``     Name of the emission line                                                      
RESTWAVE   ``float64``  Rest wavelength of the emission line                                           
ACTION     ``str_``     Action to take for this emission line; see :ref:`emission-line-modeling-action`
FLUXRATIO  ``float64``  Fixed flux ratio compared to reference line; see MODE                          
MODE       ``str_``     Modeling mode to adopt for this line; see :ref:`emission-line-modeling-mode`   
PROFILE    ``str_``     Name of the parameterization used for the instrinsic line profile              
NCOMP      ``int64``    Number of components; never used!                                              
=========  ===========  ===============================================================================

