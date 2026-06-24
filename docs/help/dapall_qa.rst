.. code-block:: console

    $ dapall_qa -h
    usage: dapall_qa [-h] [--drpver DRPVER] [--dapver DAPVER]
                     [--redux_path REDUX_PATH] [--analysis_path ANALYSIS_PATH]
                     [--plan_file PLAN_FILE] [--normal_backend]
    
    Construct QA plots using the DAPall metadata
    
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
    