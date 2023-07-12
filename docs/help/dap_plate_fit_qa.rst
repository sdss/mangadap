.. code-block:: console

    $ dap_plate_fit_qa -h
    usage: dap_plate_fit_qa [-h] [--drpver DRPVER] [--dapver DAPVER]
                            [--redux_path REDUX_PATH]
                            [--analysis_path ANALYSIS_PATH] [--plan_file PLAN_FILE]
                            [--normal_backend]
                            plate
    
    Construct per-plate QA plots
    
    positional arguments:
      plate                 plate ID to process
    
    options:
      -h, --help            show this help message and exit
      --drpver DRPVER       DRP version (default: None)
      --dapver DAPVER       DAP version (default: None)
      --redux_path REDUX_PATH
                            main DRP output path (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to use
                            instead of the default (default: None)
      --normal_backend
    