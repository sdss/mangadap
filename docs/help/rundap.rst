.. code-block:: console

    $ rundap -h
    usage: rundap [-h] [--overwrite] [-v] [--quiet] [--print_version]
                  [--drpver DRPVER] [--redux_path REDUX_PATH] [--dapver DAPVER]
                  [--analysis_path ANALYSIS_PATH] [--plan_file PLAN_FILE]
                  [--platelist PLATELIST] [--ifudesignlist IFUDESIGNLIST]
                  [--list_file LIST_FILE] [--combinatorics] [--sres_ext SRES_EXT]
                  [--sres_fill SRES_FILL] [--covar_ext COVAR_EXT] [--on_disk]
                  [--can_analyze] [--log] [--no_proc] [--no_plots] [--post]
                  [--post_plots] [--label LABEL] [--nodes NODES] [--cpus CPUS]
                  [--fast QOS] [--umask UMASK] [--walltime WALLTIME] [--toughness]
                  [--create] [--submit] [--progress] [--queue QUEUE]
    
    Perform analysis of integral-field data.
    
    options:
      -h, --help            show this help message and exit
      --overwrite           if all selected, will run dap for all
                            plates/ifudesigns/modes regardless of state (default:
                            False)
      -v, --verbose         Set verbosity level for manga_dap; can be omitted and
                            set up to -vv (default: 0)
      --quiet               suppress screen output (default: False)
      --print_version       print DAP version and stop (default: False)
      --drpver DRPVER       MaNGA DRP version for analysis; $MANGADRP_VER by default
                            (default: None)
      --redux_path REDUX_PATH
                            main DRP output path (default: None)
      --dapver DAPVER       optional output version, different from product version.
                            This *only* affects the output directory structure. It
                            does *not* select the version of the DAP to use.
                            (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to use
                            instead of the default (default: None)
      --platelist PLATELIST
                            set list of plates to reduce (default: None)
      --ifudesignlist IFUDESIGNLIST
                            set list of ifus to reduce (default: None)
      --list_file LIST_FILE
                            A file with the list of plates and ifudesigns to analyze
                            (default: None)
      --combinatorics       force execution of all permutations of the provided
                            lists (default: False)
      --sres_ext SRES_EXT   Spectral resolution extension to use. Default set by
                            MaNGADataCube class. (default: None)
      --sres_fill SRES_FILL
                            If present, use interpolation to fill any masked pixels
                            in the spectral resolution vectors. Default set by
                            MaNGADataCube class. (default: None)
      --covar_ext COVAR_EXT
                            Use this extension to define the spatial correlation
                            matrix. Default set by MaNGADataCube class. (default:
                            None)
      --on_disk             When using the DRPall file to collate the data for input
                            to the DAP, search for available DRP files on disk
                            instead of using the DRPall file content. (default:
                            False)
      --can_analyze         Only construct script files for datacubes that
                            can/should be analyzed by the DAP. See :func:`~mangadap.
                            survey.drpcomplete.DRPComplete.can_analyze`. (default:
                            False)
      --log                 Have the main DAP executable produce a log file
                            (default: False)
      --no_proc             Do NOT perform the main DAP processing steps (default:
                            False)
      --no_plots            Do NOT create QA plots (default: False)
      --post                Create/Submit the post-processing scripts (default:
                            False)
      --post_plots          Create/Submit the post-processing plotting scripts
                            (default: False)
      --label LABEL         label for cluster job (default: mangadap)
      --nodes NODES         number of nodes to use in cluster (default: 1)
      --cpus CPUS           number of cpus to use per node. Default is to use all
                            available; otherwise, set to minimum of provided number
                            and number of processors per node (default: None)
      --fast QOS            qos state (default: None)
      --umask UMASK         umask bit for cluster job (default: 0027)
      --walltime WALLTIME   walltime for cluster job (default: 240:00:00)
      --toughness           turn off hard keyword for cluster submission (default:
                            True)
      --create              use the pbs package to create the cluster scripts
                            (default: False)
      --submit              submit the scripts to the cluster (default: False)
      --progress            instead of closing the script, report the progress of
                            the analysis on the cluster; this is required if you
                            want to submit the DAPall script immediately after
                            completing the individual cube analysis (default: False)
      --queue QUEUE         set the destination queue (default: None)
    