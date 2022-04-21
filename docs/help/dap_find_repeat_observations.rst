.. code-block:: console

    $ dap_find_repeat_observations -h
    usage: dap_find_repeat_observations [-h] [--drpver DRPVER] [--dapver DAPVER]
                                        [--analysis_path ANALYSIS_PATH]
                                        [--dapall DAPALL]
                                        [--output_file OUTPUT_FILE] [-o]
    
    optional arguments:
      -h, --help            show this help message and exit
      --drpver DRPVER       DRP version (default: None)
      --dapver DAPVER       DAP version (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
      --dapall DAPALL       full path to DAPall file (default: None)
      --output_file OUTPUT_FILE
                            output file name, including path; no file is written
                            by default (default: None)
      -o, --overwrite       Overwrite any existing files. (default: False)
    