Class Instantiation: :class:`mangadap.util.drpfits.DRPQuality3DBitMask`

===========  ===  ===========================================================
Key          Bit  Description                                                
===========  ===  ===========================================================
RETIRED      0    Bit retired from use                                       
BADDEPTH     1    IFU does not reach target depth                            
SKYSUBBAD    2    Bad sky subtraction in one or more frames                  
HIGHSCAT     3    High scattered light in one or more frames                 
BADASTROM    4    Bad astrometry in one or more frames                       
VARIABLELSF  5    LSF varies signif. between component spectra               
BADOMEGA     6    Omega greater than threshhold in one or more sets          
BADSET       7    One or more sets are bad                                   
BADFLUX      8    Bad flux calibration                                       
BADPSF       9    PSF estimate may be bad                                    
MANYDEAD     10   Many dead fibers                                           
RETIRED2     11   Bit retired, moved into MASTAR_QUAL instead                
BLOWTORCH    12   Blowtorch artifact in one or more frames                   
SEVEREBT     13   Severe blowtorch artifact                                  
UNUSUAL      14   Cube slightly unusual in some way but good for most science
CRITICAL     30   Critical failure in one or more frames                     
===========  ===  ===========================================================

