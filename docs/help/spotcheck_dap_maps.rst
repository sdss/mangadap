.. code-block:: console

    $ spotcheck_dap_maps -h
    [1;34musage: [0m[1;35mspotcheck_dap_maps[0m [[32m-h[0m] ([32m-c [33mCONFIG[0m | [32m-f [33mCUBEFILE[0m)
                              [[36m--cube_module [33m[CUBE_MODULE ...][0m] [[32m-p [33mPLAN[0m]
                              [[36m--plan_module [33m[PLAN_MODULE ...][0m] [[32m-o [33mOUTPUT_PATH[0m]
                              [[32m-b [33mBEAM[0m] [[36m--normal_backend[0m]
    
    Construct a QA plot to spotcheck DAP results
    
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
      [1;32m-p[0m, [1;36m--plan[0m [1;33mPLAN[0m       SDSS parameter file with analysis plan. If not provided,
                            a default plan is used. (default: None)
      [1;36m--plan_module[0m [1;33m[PLAN_MODULE ...][0m
                            The name of the module used to define the analysis plan
                            and the output paths. (default:
                            mangadap.config.manga.MaNGAAnalysisPlan)
      [1;32m-o[0m, [1;36m--output_path[0m [1;33mOUTPUT_PATH[0m
                            Top-level directory for the DAP output files; default
                            path is set by the provided analysis plan object (see
                            plan_module). (default: None)
      [1;32m-b[0m, [1;36m--beam[0m [1;33mBEAM[0m       Beam FWHM for plot. (default: None)
      [1;36m--normal_backend[0m
    