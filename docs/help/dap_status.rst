.. code-block:: console

    $ dap_status -h
    usage: dap_status [-h] [--logdir LOGDIR] [--drpver DRPVER] [--dapver DAPVER]
                      [--analysis_path ANALYSIS_PATH]
                      plan_file
    
    positional arguments:
      plan_file             parameter file with the MaNGA DAP execution plan to
                            use instead of the default
    
    optional arguments:
      -h, --help            show this help message and exit
      --logdir LOGDIR       Log path (e.g., 31Jan2019T19.16.28UTC) (default: None)
      --drpver DRPVER       DRP version (default: None)
      --dapver DAPVER       DAP version (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
    