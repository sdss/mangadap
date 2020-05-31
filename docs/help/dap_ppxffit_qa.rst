.. code-block:: console

    $ dap_ppxffit_qa -h
    usage: dap_ppxffit_qa [-h] [--drpver DRPVER] [--redux_path REDUX_PATH]
                          [--dapver DAPVER] [--dap_src DAP_SRC]
                          [--analysis_path ANALYSIS_PATH] [--plan_file PLAN_FILE]
                          [--normal_backend]
                          [--template_flux_file TEMPLATE_FLUX_FILE]
                          plate ifudesign
    
    positional arguments:
      plate                 plate ID to process
      ifudesign             IFU design to process
    
    optional arguments:
      -h, --help            show this help message and exit
      --drpver DRPVER       DRP version (default: None)
      --redux_path REDUX_PATH
                            main DRP output path (default: None)
      --dapver DAPVER       DAP version (default: None)
      --dap_src DAP_SRC     Top-level directory with the DAP source code; defaults
                            to $MANGADAP_DIR (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to
                            use instead of the default (default: None)
      --normal_backend
      --template_flux_file TEMPLATE_FLUX_FILE
                            template renormalization flux file. Will attempt to
                            read default if not provided. If no file is provided
                            and the default file does not exist, no
                            renormalization of the templates is performed.
                            (default: None)
    