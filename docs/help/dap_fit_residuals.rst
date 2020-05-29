.. code-block:: console

    $ dap_fit_residuals -h
    usage: dap_fit_residuals [-h] [--drpver DRPVER] [--dapver DAPVER]
                             [--dap_src DAP_SRC] [--redux_path REDUX_PATH]
                             [--analysis_path ANALYSIS_PATH]
                             [--plan_file PLAN_FILE] [--daptype DAPTYPE]
                             [--normal_backend]
                             plate ifudesign
    
    positional arguments:
      plate                 plate ID to process
      ifudesign             IFU design to process
    
    optional arguments:
      -h, --help            show this help message and exit
      --drpver DRPVER       DRP version (default: None)
      --dapver DAPVER       DAP version (default: None)
      --dap_src DAP_SRC     Top-level directory with the DAP source code; defaults
                            to $MANGADAP_DIR (default: None)
      --redux_path REDUX_PATH
                            main DRP output path (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to
                            use instead of the default (default: None)
      --daptype DAPTYPE     DAP processing type (default: None)
      --normal_backend
    