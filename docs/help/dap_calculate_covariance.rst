.. code-block:: console

    $ dap_calculate_covariance -h
    usage: dap_calculate_covariance [-h] [-n NUMCHANNELS | -w WAVELENGTH]
                                    [-d DIRECTORY_PATH]
                                    plate ifudesign output_file
    
    positional arguments:
      plate                 plate ID to process
      ifudesign             IFU design to process
      output_file           Name for output file
    
    optional arguments:
      -h, --help            show this help message and exit
      -n NUMCHANNELS, --numchannels NUMCHANNELS
                            Number of channels spread across the wavelength range
                            for which to compute the covariance matrix. A value of
                            0 forces construction of the full covariance cube. The
                            default is to calculate the covariance matrix for a
                            single channel at the central wavelength (default: 1)
      -w WAVELENGTH, --wavelength WAVELENGTH
                            Wavelength at which to compute a single covariance
                            matrix; default is the central wavelength (default:
                            None)
      -d DIRECTORY_PATH, --directory_path DIRECTORY_PATH
                            Directory with the DRP produced RSS file; default uses
                            environmental variables to define the default MaNGA
                            DRP redux path (default: None)
    