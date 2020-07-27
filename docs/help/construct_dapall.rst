.. code-block:: console

    $ construct_dapall -h
    usage: construct_dapall [-h] [--plan_file PLAN_FILE] [--drpver DRPVER]
                            [-r REDUX_PATH] [--dapver DAPVER] [-a ANALYSIS_PATH]
                            [-m METHODS [METHODS ...]] [-v] [--quiet] [--double]
    
    optional arguments:
      -h, --help            show this help message and exit
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to
                            use instead of the default (default: None)
      --drpver DRPVER       DRP version (default: None)
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Top-level directory with the DRP products; defaults to
                            $MANGA_SPECTRO_REDUX/$MANGADRP_VER (default: None)
      --dapver DAPVER       DAP version (default: None)
      -a ANALYSIS_PATH, --analysis_path ANALYSIS_PATH
                            Top-level output directory for the DAP results;
                            defaults to
                            $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
                            (default: None)
      -m METHODS [METHODS ...], --methods METHODS [METHODS ...]
                            Only include output from this DAP method designation
                            in the output (default: None)
      -v, --verbose         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      --quiet               suppress all terminal output (default: False)
      --double              Output the floating-point data in double precision
                            (default is single precision) (default: True)
    