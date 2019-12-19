
=================  =============  =======  =======  =====================================================================================================================================================================
Key                Type           Options  Default  Description                                                                                                                                                          
=================  =============  =======  =======  =====================================================================================================================================================================
``key``            str            ..       ..       Keyword to distinguish the assessment method.                                                                                                                        
``waverange``      ndarray, list  ..       ..       A two-element vector with the starting and ending wavelength (angstroms in **vacuum**) within which to calculate the signal-to-noise                                 
``response_func``  ndarray, list  ..       ..       A two-column array with a response function to use for the S/N calculation.  The columns must br the wavelength and amplitude of the response function, respectively.
``covariance``     bool           ..       ..       Type of covariance measurement to produce.                                                                                                                           
=================  =============  =======  =======  =====================================================================================================================================================================

