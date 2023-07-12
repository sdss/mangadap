.. code-block:: console

    $ manga_dap -h
    usage: manga_dap [-h] (-c CONFIG | -f CUBEFILE)
                     [--cube_module [CUBE_MODULE ...]] [-p PLAN]
                     [--plan_module [PLAN_MODULE ...]] [--dbg] [--log LOG] [-v]
                     [-o OUTPUT_PATH]
    
    Perform analysis of integral-field data.
    
    options:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Configuration file used to instantiate the relevant
                            DataCube derived class. (default: None)
      -f CUBEFILE, --cubefile CUBEFILE
                            Name of the file with the datacube data. Must be
                            possible to instantiate the relevant DataCube derived
                            class directly from the file only. (default: None)
      --cube_module [CUBE_MODULE ...]
                            The name of the module that contains the DataCube
                            derived class used to read the data. (default:
                            mangadap.datacube.MaNGADataCube)
      -p PLAN, --plan PLAN  TOML file with analysis plan. If not provided, a default
                            plan is used. (default: None)
      --plan_module [PLAN_MODULE ...]
                            The name of the module used to define the analysis plan
                            and the output paths. (default:
                            mangadap.config.manga.MaNGAAnalysisPlan)
      --dbg                 Run manga_dap in debug mode (default: False)
      --log LOG             File name for runtime log (default: None)
      -v, --verbose         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      -o OUTPUT_PATH, --output_path OUTPUT_PATH
                            Top-level directory for the DAP output files; default
                            path is set by the provided analysis plan object (see
                            plan_module). (default: None)
    