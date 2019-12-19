
=======================  ==========  =======  =======  =================================================================================================
Key                      Type        Options  Default  Description                                                                                      
=======================  ==========  =======  =======  =================================================================================================
``key``                  str         ..       ..       Keyword used to distinguish between different spectral-index databases.                          
``minimum_snr``          int, float  ..       ..       Minimum S/N of spectrum to fit                                                                   
``fwhm``                 int, float  ..       ..       Resolution FWHM in angstroms at which to make the measurements.                                  
``compute_corrections``  bool        ..       ..       Flag to compute dispersion corrections to indices.  Dispersion corrections are always calculated!
``artifacts``            str         ..       ..       String identifying the artifact database to use                                                  
``absindex``             str         ..       ..       String identifying the absorption-index database to use                                          
``bandhead``             str         ..       ..       String identifying the bandhead-index database to use                                            
=======================  ==========  =======  =======  =================================================================================================

