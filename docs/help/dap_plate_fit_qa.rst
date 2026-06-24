.. code-block:: console

    $ dap_plate_fit_qa -h
    [1;34musage: [0m[1;35mdap_plate_fit_qa[0m [[32m-h[0m] [[36m--drpver [33mDRPVER[0m] [[36m--dapver [33mDAPVER[0m]
                            [[36m--redux_path [33mREDUX_PATH[0m]
                            [[36m--analysis_path [33mANALYSIS_PATH[0m] [[36m--plan_file [33mPLAN_FILE[0m]
                            [[36m--normal_backend[0m]
                            [32mplate[0m
    
    Construct per-plate QA plots
    
    [1;34mpositional arguments:[0m
      [1;32mplate[0m                 plate ID to process
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;36m--drpver[0m [1;33mDRPVER[0m       DRP version (default: None)
      [1;36m--dapver[0m [1;33mDAPVER[0m       DAP version (default: None)
      [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            main DRP output path (default: None)
      [1;36m--analysis_path[0m [1;33mANALYSIS_PATH[0m
                            main DAP output path (default: None)
      [1;36m--plan_file[0m [1;33mPLAN_FILE[0m
                            parameter file with the MaNGA DAP execution plan to use
                            instead of the default (default: None)
      [1;36m--normal_backend[0m
    