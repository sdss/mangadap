.. code-block:: console

    $ construct_dapall -h
    usage: construct_dapall [-h] [--drpver DRPVER] [-r REDUX_PATH] [--dapver DAPVER]
                            [-a ANALYSIS_PATH] [--plan_file PLAN_FILE]
                            [-m [METHODS ...]] [-v] [--quiet] [--double]
    
    Compile metadata for the DAPall file
    
    options:
      -h, --help            show this help message and exit
      --drpver DRPVER       DRP version (default: None)
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Top-level directory with the DRP products; defaults to
                            $MANGA_SPECTRO_REDUX/$MANGADRP_VER (default: None)
      --dapver DAPVER       DAP version (default: None)
      -a ANALYSIS_PATH, --analysis_path ANALYSIS_PATH
                            Top-level output directory for the DAP results; defaults
                            to $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
                            (default: None)
      --plan_file PLAN_FILE
                            toml file with the MaNGA DAP execution plan to use
                            instead of the default (default: None)
      -m [METHODS ...], --methods [METHODS ...]
                            Only include output for these DAP method designations in
                            the output (default: None)
      -v, --verbose         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      --quiet               suppress all terminal output (default: False)
      --double              Output the floating-point data in double precision
                            (default is single precision) (default: True)
    