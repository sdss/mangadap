
======================  ============  =======  ============  ===============================================================================================================================================================
Key                     Type          Options  Default       Description                                                                                                                                                    
======================  ============  =======  ============  ===============================================================================================================================================================
``key``                 str           ..       ``SPX``       Keyword used to distinguish between different spatial binning schemes.                                                                                         
``galactic_reddening``  str           ..       ``ODonnell``  The string identifier for the Galactic extinction curve to use.  See :func:`~mangadap.util.extinction.GalacticExtinction.valid_forms` for the available curves.
``galactic_rv``         int, float    ..       3.1           Ratio of V-band extinction to the B-V reddening.                                                                                                               
``minimum_snr``         int, float    ..       1.0           Minimum S/N of spectra to include in any bin.                                                                                                                  
``minimum_frac``        int, float    ..       0.8           Minimum fraction of unmasked pixels in each spectrum included in any bin.                                                                                      
``binpar``              ParSet, dict  ..       ..            The spatial-binning parameters.                                                                                                                                
``binclass``            Undefined     ..       ..            Instance of the spatial-binning class.  Needed in case binfunc is a non-static member function of the class.                                                   
``binfunc``             Undefined     ..       ..            The spatial-binning function that determines which spectra go into each bin.                                                                                   
``stackpar``            ParSet, dict  ..       ..            The spectral-stacking parameter set.                                                                                                                           
``stackclass``          Undefined     ..       ..            Instance of spectral-stacking class to use.  Needed in case stackfunc is a non-static member function of the class.                                            
``stackfunc``           Undefined     ..       ..            The spectral-stacking function that stacks the spectra in a given bin.                                                                                         
``overwrite``           bool          ..       False         If the output file already exists, redo all the calculations and overwrite it.                                                                                 
======================  ============  =======  ============  ===============================================================================================================================================================

