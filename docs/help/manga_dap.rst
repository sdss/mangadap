.. code-block:: console

    $ manga_dap -h
    usage: manga_dap [-h] (-c CONFIG | -f CUBEFILE) [-p PLAN] [-m CUBE_MODULE]
                     [-o CUBE_OBJECT] [--dbg] [--log LOG] [-v] [--drpver DRPVER]
                     [-r REDUX_PATH] [-d DIRECTORY_PATH] [--dapver DAPVER]
                     [-a ANALYSIS_PATH]
    
    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Configuration file used to instantiate the relevant
                            DataCube derived class. (default: None)
      -f CUBEFILE, --cubefile CUBEFILE
                            Name of the file with the datacube data. Must be
                            possible to instantiate the relevant DataCube derived
                            class directly from the file only. (default: None)
      -p PLAN, --plan PLAN  SDSS parameter file with analysis plan. If not
                            provided, a default plan is used. (default: None)
      -m CUBE_MODULE, --cube_module CUBE_MODULE
                            The name of the module that contains the DataCube
                            derived class. (default: mangadap.datacube)
      -o CUBE_OBJECT, --cube_object CUBE_OBJECT
                            The name of the DataCube derived class object.
                            (default: MaNGADataCube)
      --dbg                 Run manga_dap in debug mode (default: False)
      --log LOG             File name for runtime log (default: None)
      -v, --verbose         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      --drpver DRPVER       DRP version (default: None)
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Top-level directory with the DRP products; defaults to
                            $MANGA_SPECTRO_REDUX/$MANGADRP_VER (default: None)
      -d DIRECTORY_PATH, --directory_path DIRECTORY_PATH
                            Path directly to directory with DRP file to analyze
                            (default: None)
      --dapver DAPVER       DAP version (default: None)
      -a ANALYSIS_PATH, --analysis_path ANALYSIS_PATH
                            Top-level output directory for the DAP results;
                            defaults to
                            $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
                            (default: None)
    