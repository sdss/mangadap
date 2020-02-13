
======================  ============  ====================  =======  ========================================================================================================================================================================================
Key                     Type          Options               Default  Description                                                                                                                                                                             
======================  ============  ====================  =======  ========================================================================================================================================================================================
``key``                 str           ..                    ..       Keyword used to distinguish between different spatial binning schemes.                                                                                                                  
``galactic_reddening``  str           ..                    ..       The string identifier for the Galactic extinction curve to use.  See :func:`mangadap.util.extinction.GalacticExtinction.valid_forms` for the available curves.  Default is ``ODonnell``.
``galactic_rv``         int, float    ..                    ..       Ratio of V-band extinction to the B-V reddening.  Default is 3.1.                                                                                                                       
``minimum_snr``         int, float    ..                    ..       Minimum S/N of spectra to include in any bin.                                                                                                                                           
``binpar``              ParSet, dict  ..                    ..       The parameter set defining how to place each spectrum in a bin.                                                                                                                         
``binclass``            Undefined     ..                    ..       Instance of class object to use for the binning.  Needed in case binfunc is a non-static member function of the class.                                                                  
``binfunc``             Undefined     ..                    ..       The function that determines which spectra go into each bin.                                                                                                                            
``stackpar``            ParSet, dict  ..                    ..       The parameter set defining how to stack the spectra in each bin.                                                                                                                        
``stackclass``          Undefined     ..                    ..       Instance of class object to used to stack the spectra.  Needed in case stackfunc is a non-static member function of the class.                                                          
``stackfunc``           Undefined     ..                    ..       The function that stacks the spectra in a given bin.                                                                                                                                    
``spec_res``            str           ``spaxel``, ``cube``  ..       Keyword defining the treatment of the spectral resolution.  See :func:`SpatiallyBinnedSpectra.spectral_resolution_options` for a list of the options.                                   
``prepixel_sres``       bool          ..                    ..       Use the prepixelized version of the LSF measurements.                                                                                                                                   
======================  ============  ====================  =======  ========================================================================================================================================================================================
