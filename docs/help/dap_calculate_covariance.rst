.. code-block:: console

    $ dap_calculate_covariance -h
    [1;34musage: [0m[1;35mdap_calculate_covariance[0m [[32m-h[0m] [[32m-n [33mNUMCHANNELS[0m | [32m-w [33mWAVELENGTH[0m]
                                    [[32m-d [33mDIRECTORY_PATH[0m]
                                    [32mplate[0m [32mifudesign[0m [32moutput_file[0m
    
    Calculate covariance in a MaNGA datacube
    
    [1;34mpositional arguments:[0m
      [1;32mplate[0m                 plate ID to process
      [1;32mifudesign[0m             IFU design to process
      [1;32moutput_file[0m           Name for output file
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-n[0m, [1;36m--numchannels[0m [1;33mNUMCHANNELS[0m
                            Number of channels spread across the wavelength range
                            for which to compute the covariance matrix. A value of 0
                            forces construction of the full covariance cube. The
                            default is to calculate the covariance matrix for a
                            single channel at the central wavelength (default: 1)
      [1;32m-w[0m, [1;36m--wavelength[0m [1;33mWAVELENGTH[0m
                            Wavelength at which to compute a single covariance
                            matrix; default is the central wavelength (default:
                            None)
      [1;32m-d[0m, [1;36m--directory_path[0m [1;33mDIRECTORY_PATH[0m
                            Directory with the DRP produced RSS file; default uses
                            environmental variables to define the default MaNGA DRP
                            redux path (default: None)
    