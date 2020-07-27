.. code-block:: console

    $ write_dap_config -h
    usage: write_dap_config [-h] (-c DRPCOMPLETE | -a DRPALL)
                            [--sres_ext SRES_EXT] [--sres_fill SRES_FILL]
                            [--covar_ext COVAR_EXT] [--drpver DRPVER]
                            [--redux_path REDUX_PATH]
                            [--directory_path DIRECTORY_PATH] [-o]
                            plate ifudesign ofile
    
    positional arguments:
      plate                 Plate number
      ifudesign             IFU design number
      ofile                 Output file name
    
    optional arguments:
      -h, --help            show this help message and exit
      -c DRPCOMPLETE, --drpcomplete DRPCOMPLETE
                            DRP complete fits file (default: None)
      -a DRPALL, --drpall DRPALL
                            DRPall fits file (default: None)
      --sres_ext SRES_EXT   Spectral resolution extension to use. Default set by
                            MaNGADataCube class. (default: None)
      --sres_fill SRES_FILL
                            If present, use interpolation to fill any masked
                            pixels in the spectral resolution vectors. Default set
                            by MaNGADataCube class. (default: None)
      --covar_ext COVAR_EXT
                            Use this extension to define the spatial correlation
                            matrix. Default set by MaNGADataCube class. (default:
                            None)
      --drpver DRPVER       DRP version. Default set by MaNGADataCube class.
                            (default: None)
      --redux_path REDUX_PATH
                            Path to the top-level DRP reduction directory. Default
                            set by MaNGADataCube class. (default: None)
      --directory_path DIRECTORY_PATH
                            Exact path to the directory with the MaNGA DRP
                            datacube. The name of the file itself must match the
                            nominal MaNGA DRP naming convention. Default set by
                            MaNGADataCube class. (default: None)
      -o, --overwrite       Overwrite any existing files. (default: False)
    