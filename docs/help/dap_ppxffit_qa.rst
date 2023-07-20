.. code-block:: console

    $ dap_ppxffit_qa -h
    usage: dap_ppxffit_qa [-h] (-c CONFIG | -f CUBEFILE)
                          [--cube_module [CUBE_MODULE ...]] [-p PLAN]
                          [--plan_module [PLAN_MODULE ...]] [-o OUTPUT_PATH]
                          [-b BEAM] [--normal_backend]
                          [--template_flux_file TEMPLATE_FLUX_FILE]
    
    Construct QA plots for pPXF fit.
    
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
      -p PLAN, --plan PLAN  SDSS parameter file with analysis plan. If not provided,
                            a default plan is used. (default: None)
      --plan_module [PLAN_MODULE ...]
                            The name of the module used to define the analysis plan
                            and the output paths. (default:
                            mangadap.config.manga.MaNGAAnalysisPlan)
      -o OUTPUT_PATH, --output_path OUTPUT_PATH
                            Top-level directory for the DAP output files; default
                            path is set by the provided analysis plan object (see
                            plan_module). (default: None)
      -b BEAM, --beam BEAM  Beam FWHM for plot. (default: None)
      --normal_backend
      --template_flux_file TEMPLATE_FLUX_FILE
                            template renormalization flux file. Will attempt to read
                            default if not provided. If no file is provided and the
                            default file does not exist, no renormalization of the
                            templates is performed. (default: None)
    