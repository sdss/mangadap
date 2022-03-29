
=================  ============  =======  ==================  =========================================================================================================
Key                Type          Options  Default             Description                                                                                              
=================  ============  =======  ==================  =========================================================================================================
``key``            str           ..       ``EFITMPL11SSPDB``  Keyword used to distinguish between different emission-line moment databases.                            
``minimum_snr``    int, float    ..       0.0                 Minimum S/N of spectrum to fit                                                                           
``mom_vel_name``   str           ..       ``Ha-6564``         Name of the emission-line moments band used to set the initial velocity guess for each spaxel.           
``mom_disp_name``  str           ..       ..                  Name of the emission-line moments band used to set the initial velocity dispersion guess for each spaxel.
``fitpar``         ParSet, dict  ..       ..                  Fitting function parameters                                                                              
``fitclass``       Undefined     ..       ..                  Class used to perform the fit.                                                                           
``fitfunc``        Undefined     ..       ..                  Function or method that performs the fit; **must** be callable.                                          
``overwrite``      bool          ..       False               If the output file already exists, redo all the calculations and overwrite it.                           
=================  ============  =======  ==================  =========================================================================================================

