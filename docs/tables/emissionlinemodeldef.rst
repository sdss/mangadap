
=====================  ============  =======  =======  ======================================================================================================================================================================
Key                    Type          Options  Default  Description                                                                                                                                                           
=====================  ============  =======  =======  ======================================================================================================================================================================
``key``                str           ..       ..       Keyword used to distinguish between different emission-line moment databases.                                                                                         
``minimum_snr``        int, float    ..       ..       Minimum S/N of spectrum to fit                                                                                                                                        
``deconstruct_bins``   str           ..       ..       Method to use for deconstructing binned spectra into individual spaxels for emission-line fitting.  See :func:`~mangadap.proc.sasuke.Sasuke.deconstruct_bins_options`.
``mom_vel_name``       str           ..       ..       Name of the emission-line moments band used to set the initial velocity guess for each spaxel.                                                                        
``mom_disp_name``      str           ..       ..       Name of the emission-line moments band used to set the initial velocity dispersion guess for each spaxel.                                                             
``waverange``          list          ..       ..       Limited wavelength range to use during the fit                                                                                                                        
``artifacts``          str           ..       ..       String identifying the artifact database to use                                                                                                                       
``ism_mask``           str           ..       ..       String identifying an emission-line database used only for **masking** lines during the fit.                                                                          
``emission_lines``     str           ..       ..       String identifying the emission-line database to use                                                                                                                  
``continuum_tpl_key``  str           ..       ..       String identifying the continuum templates to use                                                                                                                     
``fitpar``             ParSet, dict  ..       ..       Fitting function parameters                                                                                                                                           
``fitclass``           Undefined     ..       ..       Class used to perform the fit.                                                                                                                                        
``fitfunc``            Undefined     ..       ..       Function or method that performs the fit; **must** be callable.                                                                                                       
=====================  ============  =======  =======  ======================================================================================================================================================================

