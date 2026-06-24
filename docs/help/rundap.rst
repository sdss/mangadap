.. code-block:: console

    $ rundap -h
    [1;34musage: [0m[1;35mrundap[0m [[32m-h[0m] [[36m--overwrite[0m] [[32m-v[0m] [[36m--quiet[0m] [[36m--print_version[0m]
                  [[36m--drpver [33mDRPVER[0m] [[36m--redux_path [33mREDUX_PATH[0m] [[36m--dapver [33mDAPVER[0m]
                  [[36m--analysis_path [33mANALYSIS_PATH[0m] [[36m--plan_file [33mPLAN_FILE[0m]
                  [[36m--platelist [33mPLATELIST[0m] [[36m--ifudesignlist [33mIFUDESIGNLIST[0m]
                  [[36m--list_file [33mLIST_FILE[0m] [[36m--combinatorics[0m] [[36m--sres_ext [33mSRES_EXT[0m]
                  [[36m--sres_fill [33mSRES_FILL[0m] [[36m--covar_ext [33mCOVAR_EXT[0m] [[36m--on_disk[0m]
                  [[36m--can_analyze[0m] [[36m--log[0m] [[36m--no_proc[0m] [[36m--no_plots[0m] [[36m--post[0m]
                  [[36m--post_plots[0m] [[36m--label [33mLABEL[0m] [[36m--nodes [33mNODES[0m] [[36m--cpus [33mCPUS[0m]
                  [[36m--fast [33mQOS[0m] [[36m--umask [33mUMASK[0m] [[36m--walltime [33mWALLTIME[0m] [[36m--toughness[0m]
                  [[36m--create[0m] [[36m--submit[0m] [[36m--progress[0m] [[36m--queue [33mQUEUE[0m]
    
    Perform analysis of integral-field data.
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;36m--overwrite[0m           if all selected, will run dap for all
                            plates/ifudesigns/modes regardless of state (default:
                            False)
      [1;32m-v[0m, [1;36m--verbose[0m         Set verbosity level for manga_dap; can be omitted and
                            set up to -vv (default: 0)
      [1;36m--quiet[0m               suppress screen output (default: False)
      [1;36m--print_version[0m       print DAP version and stop (default: False)
      [1;36m--drpver[0m [1;33mDRPVER[0m       MaNGA DRP version for analysis; $MANGADRP_VER by default
                            (default: None)
      [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            main DRP output path (default: None)
      [1;36m--dapver[0m [1;33mDAPVER[0m       optional output version, different from product version.
                            This *only* affects the output directory structure. It
                            does *not* select the version of the DAP to use.
                            (default: None)
      [1;36m--analysis_path[0m [1;33mANALYSIS_PATH[0m
                            main DAP output path (default: None)
      [1;36m--plan_file[0m [1;33mPLAN_FILE[0m
                            parameter file with the MaNGA DAP execution plan to use
                            instead of the default (default: None)
      [1;36m--platelist[0m [1;33mPLATELIST[0m
                            set list of plates to reduce (default: None)
      [1;36m--ifudesignlist[0m [1;33mIFUDESIGNLIST[0m
                            set list of ifus to reduce (default: None)
      [1;36m--list_file[0m [1;33mLIST_FILE[0m
                            A file with the list of plates and ifudesigns to analyze
                            (default: None)
      [1;36m--combinatorics[0m       force execution of all permutations of the provided
                            lists (default: False)
      [1;36m--sres_ext[0m [1;33mSRES_EXT[0m   Spectral resolution extension to use. Default set by
                            MaNGADataCube class. (default: None)
      [1;36m--sres_fill[0m [1;33mSRES_FILL[0m
                            If present, use interpolation to fill any masked pixels
                            in the spectral resolution vectors. Default set by
                            MaNGADataCube class. (default: None)
      [1;36m--covar_ext[0m [1;33mCOVAR_EXT[0m
                            Use this extension to define the spatial correlation
                            matrix. Default set by MaNGADataCube class. (default:
                            None)
      [1;36m--on_disk[0m             When using the DRPall file to collate the data for input
                            to the DAP, search for available DRP files on disk
                            instead of using the DRPall file content. (default:
                            False)
      [1;36m--can_analyze[0m         Only construct script files for datacubes that
                            can/should be analyzed by the DAP. See :func:`~mangadap.
                            survey.drpcomplete.DRPComplete.can_analyze`. (default:
                            False)
      [1;36m--log[0m                 Have the main DAP executable produce a log file
                            (default: False)
      [1;36m--no_proc[0m             Do NOT perform the main DAP processing steps (default:
                            False)
      [1;36m--no_plots[0m            Do NOT create QA plots (default: False)
      [1;36m--post[0m                Create/Submit the post-processing scripts (default:
                            False)
      [1;36m--post_plots[0m          Create/Submit the post-processing plotting scripts
                            (default: False)
      [1;36m--label[0m [1;33mLABEL[0m         label for cluster job (default: mangadap)
      [1;36m--nodes[0m [1;33mNODES[0m         number of nodes to use in cluster (default: 1)
      [1;36m--cpus[0m [1;33mCPUS[0m           number of cpus to use per node. Default is to use all
                            available; otherwise, set to minimum of provided number
                            and number of processors per node (default: None)
      [1;36m--fast[0m [1;33mQOS[0m            qos state (default: None)
      [1;36m--umask[0m [1;33mUMASK[0m         umask bit for cluster job (default: 0027)
      [1;36m--walltime[0m [1;33mWALLTIME[0m   walltime for cluster job (default: 240:00:00)
      [1;36m--toughness[0m           turn off hard keyword for cluster submission (default:
                            True)
      [1;36m--create[0m              use the pbs package to create the cluster scripts
                            (default: False)
      [1;36m--submit[0m              submit the scripts to the cluster (default: False)
      [1;36m--progress[0m            instead of closing the script, report the progress of
                            the analysis on the cluster; this is required if you
                            want to submit the DAPall script immediately after
                            completing the individual cube analysis (default: False)
      [1;36m--queue[0m [1;33mQUEUE[0m         set the destination queue (default: None)
    