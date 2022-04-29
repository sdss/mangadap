
=====================  =================  =======  =============  =============================================================================================================================================
Key                    Type               Options  Default        Description                                                                                                                                  
=====================  =================  =======  =============  =============================================================================================================================================
``key``                str                ..       ``EMOMMPL11``  Keyword used to distinguish between different emission-line moment databases.                                                                
``minimum_snr``        int, float         ..       0.0            Minimum S/N of spectrum to analyze                                                                                                           
``pixelmask``          SpectralPixelMask  ..       ..             Object used to mask spectral pixels                                                                                                          
``passbands``          str                ..       ``ELBMPL9``    Either a string identifying the emission-line bandpass filter database to use, or the direct path to the parameter file defining the database
``redo_postmodeling``  bool               ..       True           Redo the moment measurements after the emission-line modeling has been performed                                                             
``fit_vel_name``       str                ..       ``Ha-6564``    The name of the emission line used to set the redshift of each spaxel used to set the observed wavelength of the bandpasses.                 
``overwrite``          bool               ..       False          If the output file already exists, redo all the calculations and overwrite it.                                                               
=====================  =================  =======  =============  =============================================================================================================================================

