
====================  ==========  =======  =======  =========================================================================================================
Key                   Type        Options  Default  Description                                                                                              
====================  ==========  =======  =======  =========================================================================================================
``key``               str         ..       ..       Keyword to distinguish the template library.                                                             
``file_search``       str         ..       ..       Search string used by glob to find the 1D fits spectra to include in the template library.               
``fwhm``              int, float  ..       ..       FWHM of the resolution element in angstroms.                                                             
``sres_ext``          str         ..       ..       Extension in the fits files with measurements of the spectral resolution as a function of wavelength.    
``in_vacuum``         bool        ..       ..       Flag that the wavelengths of the spectra are in vacuum, not air.                                         
``wave_limit``        ndarray     ..       ..       Two-element array with the starting and ending wavelengths for the valid spectral range of the templates.
``lower_flux_limit``  int, float  ..       ..       Minimum valid flux in the template spectra.                                                              
``log10``             bool        ..       ..       Flag that the template spectra have been binned logarithmically in wavelength.                           
====================  ==========  =======  =======  =========================================================================================================

