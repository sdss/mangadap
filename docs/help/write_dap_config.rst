.. code-block:: console

    $ write_dap_config -h
    [1;34musage: [0m[1;35mwrite_dap_config[0m [[32m-h[0m] ([32m-c [33mDRPCOMPLETE[0m | [32m-a [33mDRPALL[0m) [[36m--sres_ext [33mSRES_EXT[0m]
                            [[36m--sres_fill [33mSRES_FILL[0m] [[36m--covar_ext [33mCOVAR_EXT[0m]
                            [[36m--drpver [33mDRPVER[0m] [[36m--redux_path [33mREDUX_PATH[0m]
                            [[36m--directory_path [33mDIRECTORY_PATH[0m] [[32m-o[0m]
                            [32mplate[0m [32mifudesign[0m [32mofile[0m
    
    Generate a DAP input configuration file
    
    [1;34mpositional arguments:[0m
      [1;32mplate[0m                 Plate number
      [1;32mifudesign[0m             IFU design number
      [1;32mofile[0m                 Output file name
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-c[0m, [1;36m--drpcomplete[0m [1;33mDRPCOMPLETE[0m
                            DRP complete fits file (default: None)
      [1;32m-a[0m, [1;36m--drpall[0m [1;33mDRPALL[0m   DRPall fits file (default: None)
      [1;36m--sres_ext[0m [1;33mSRES_EXT[0m   Spectral resolution extension to use. Default set by
                            MaNGADataCube class. (default: None)
      [1;36m--sres_fill[0m [1;33mSRES_FILL[0m
                            If present, use interpolation to fill any masked pixels
                            in the spectral resolution vectors. Default set by
                            MaNGADataCube class. (default: None)
      [1;36m--covar_ext[0m [1;33mCOVAR_EXT[0m
                            Use this extension to define the spatial correlation
                            matrix. Default set by MaNGADataCube class. (default:
                            None)
      [1;36m--drpver[0m [1;33mDRPVER[0m       DRP version. Default set by MaNGADataCube class.
                            (default: None)
      [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            Path to the top-level DRP reduction directory. Default
                            set by MaNGADataCube class. (default: None)
      [1;36m--directory_path[0m [1;33mDIRECTORY_PATH[0m
                            Exact path to the directory with the MaNGA DRP datacube.
                            The name of the file itself must match the nominal MaNGA
                            DRP naming convention. Default set by MaNGADataCube
                            class. (default: None)
      [1;32m-o[0m, [1;36m--overwrite[0m       Overwrite any existing files. (default: False)
    