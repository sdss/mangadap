Class Instantiation: :class:`mangadap.dapfits.DAPMapsBitMask`

============  ===  ==============================================================================
Key           Bit  Description                                                                   
============  ===  ==============================================================================
NOCOV         0    No coverage in this spaxel (propagated from DRP)                              
LOWCOV        1    Low coverage in this spaxel (propagated from DRP)                             
DEADFIBER     2    Major contributing fiber is dead (propagated from DRP)                        
FORESTAR      3    Foreground star (propagated from DRP)                                         
NOVALUE       4    Spaxel was not included in analysis because it did not meet selection criteria
UNRELIABLE    5    Value is deemed unreliable                                                    
MATHERROR     6    Mathematical error in computing value                                         
FITFAILED     7    Fit optimization algorithm failed                                             
NEARBOUND     8    Fitted parameter is at or too near an imposed boundary                        
NOCORRECTION  9    Appropriate correction for this parameter is not available                    
MULTICOMP     10   Multi-component velocity features present                                     
DONOTUSE      30   Do not use this spaxel for science                                            
============  ===  ==============================================================================

