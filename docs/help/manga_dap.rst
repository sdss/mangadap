.. code-block:: console

    $ manga_dap -h
    [1;34musage: [0m[1;35mmanga_dap[0m [[32m-h[0m] ([32m-c [33mCONFIG[0m | [32m-f [33mCUBEFILE[0m)
                     [[36m--cube_module [33m[CUBE_MODULE ...][0m] [[32m-p [33mPLAN[0m]
                     [[36m--plan_module [33m[PLAN_MODULE ...][0m] [[36m--dbg[0m] [[36m--log [33mLOG[0m] [[32m-v[0m]
                     [[32m-o [33mOUTPUT_PATH[0m]
    
    Perform analysis of integral-field data.
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-c[0m, [1;36m--config[0m [1;33mCONFIG[0m   Configuration file used to instantiate the relevant
                            DataCube derived class. (default: None)
      [1;32m-f[0m, [1;36m--cubefile[0m [1;33mCUBEFILE[0m
                            Name of the file with the datacube data. Must be
                            possible to instantiate the relevant DataCube derived
                            class directly from the file only. (default: None)
      [1;36m--cube_module[0m [1;33m[CUBE_MODULE ...][0m
                            The name of the module that contains the DataCube
                            derived class used to read the data. (default:
                            mangadap.datacube.MaNGADataCube)
      [1;32m-p[0m, [1;36m--plan[0m [1;33mPLAN[0m       TOML file with analysis plan. If not provided, a default
                            plan is used. (default: None)
      [1;36m--plan_module[0m [1;33m[PLAN_MODULE ...][0m
                            The name of the module used to define the analysis plan
                            and the output paths. (default:
                            mangadap.config.manga.MaNGAAnalysisPlan)
      [1;36m--dbg[0m                 Run manga_dap in debug mode (default: False)
      [1;36m--log[0m [1;33mLOG[0m             File name for runtime log (default: None)
      [1;32m-v[0m, [1;36m--verbose[0m         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      [1;32m-o[0m, [1;36m--output_path[0m [1;33mOUTPUT_PATH[0m
                            Top-level directory for the DAP output files; default
                            path is set by the provided analysis plan object (see
                            plan_module). (default: None)
    