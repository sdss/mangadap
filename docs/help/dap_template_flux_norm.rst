.. code-block:: console

    $ dap_template_flux_norm -h
    usage: dap_template_flux_norm [-h] [-v VELSCALE_RATIO] [-s STEP] [-l]
                                  key ofile
    
    positional arguments:
      key                   DAP template-library key
      ofile                 Output file name
    
    optional arguments:
      -h, --help            show this help message and exit
      -v VELSCALE_RATIO, --velscale_ratio VELSCALE_RATIO
                            velocity scale ratio between the template pixel and
                            the nominal pixel width (default: 4)
      -s STEP, --step STEP  linear or log-linear spectral pixel width (default:
                            0.0001)
      -l, --linear          pixel sampling should be linear (default is log)
                            (default: False)
    