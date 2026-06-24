.. code-block:: console

    $ construct_dapall -h
    [1;34musage: [0m[1;35mconstruct_dapall[0m [[32m-h[0m] [[36m--drpver [33mDRPVER[0m] [[32m-r [33mREDUX_PATH[0m] [[36m--dapver [33mDAPVER[0m]
                            [[32m-a [33mANALYSIS_PATH[0m] [[36m--plan_file [33mPLAN_FILE[0m]
                            [[32m-m [33m[METHODS ...][0m] [[32m-v[0m] [[36m--quiet[0m] [[36m--double[0m]
    
    Compile metadata for the DAPall file
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;36m--drpver[0m [1;33mDRPVER[0m       DRP version (default: None)
      [1;32m-r[0m, [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            Top-level directory with the DRP products; defaults to
                            $MANGA_SPECTRO_REDUX/$MANGADRP_VER (default: None)
      [1;36m--dapver[0m [1;33mDAPVER[0m       DAP version (default: None)
      [1;32m-a[0m, [1;36m--analysis_path[0m [1;33mANALYSIS_PATH[0m
                            Top-level output directory for the DAP results; defaults
                            to $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
                            (default: None)
      [1;36m--plan_file[0m [1;33mPLAN_FILE[0m
                            toml file with the MaNGA DAP execution plan to use
                            instead of the default (default: None)
      [1;32m-m[0m, [1;36m--methods[0m [1;33m[METHODS ...][0m
                            Only include output for these DAP method designations in
                            the output (default: None)
      [1;32m-v[0m, [1;36m--verbose[0m         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      [1;36m--quiet[0m               suppress all terminal output (default: False)
      [1;36m--double[0m              Output the floating-point data in double precision
                            (default is single precision) (default: True)
    