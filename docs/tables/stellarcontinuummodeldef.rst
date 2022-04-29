
===============  ============  =======  ===========  =============================================================================================================================================================
Key              Type          Options  Default      Description                                                                                                                                                  
===============  ============  =======  ===========  =============================================================================================================================================================
``key``          str           ..       ``MILESHC``  Keyword used to distinguish between different spatial binning schemes.                                                                                       
``minimum_snr``  int, float    ..       1.0          Minimum S/N of spectrum to fit                                                                                                                               
``fitpar``       ParSet, dict  ..       ..           Any additional parameters, aside from the spectra themselves, required by the fitting function.  If None, uses a default pPXF fit.                           
``fitclass``     Undefined     ..       ..           Instance of class object to use for the model fitting.  Needed in case fitfunc is a non-static member function of a class.  If None, uses a default pPXF fit.
``fitfunc``      Undefined     ..       ..           The function that models the spectra.  If None, uses a default pPXF fit and anything provided for fitpar and fitclass are ignored!                           
``overwrite``    bool          ..       False        If the output file already exists, redo all the calculations and overwrite it.                                                                               
===============  ============  =======  ===========  =============================================================================================================================================================

