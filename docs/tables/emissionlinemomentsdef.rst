
=====================  ==========  =======  =======  ============================================================================================================================
Key                    Type        Options  Default  Description                                                                                                                 
=====================  ==========  =======  =======  ============================================================================================================================
``key``                str         ..       ..       Keyword used to distinguish between different emission-line moment databases.                                               
``minimum_snr``        int, float  ..       ..       Minimum S/N of spectrum to fit                                                                                              
``artifacts``          str         ..       ..       String identifying the artifact database to use                                                                             
``passbands``          str         ..       ..       String identifying the emission-line bandpass filter database to use                                                        
``redo_postmodeling``  bool        ..       ..       Redo the moment measurements after the emission-line modeling has been performed                                            
``fit_vel_name``       str         ..       ..       The name of the emission line used to set the redshift of each spaxel used to set the observed wavelength of the bandpasses.
=====================  ==========  =======  =======  ============================================================================================================================

