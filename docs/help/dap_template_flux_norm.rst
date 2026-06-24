.. code-block:: console

    $ dap_template_flux_norm -h
    [1;34musage: [0m[1;35mdap_template_flux_norm[0m [[32m-h[0m] [[32m-v [33mVELSCALE_RATIO[0m] [[32m-s [33mSTEP[0m] [[32m-l[0m] [32mkey[0m [32mofile[0m
    
    Compute flux normalizations for DAP template set
    
    [1;34mpositional arguments:[0m
      [1;32mkey[0m                   DAP template-library key
      [1;32mofile[0m                 Output file name
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--velscale_ratio[0m [1;33mVELSCALE_RATIO[0m
                            velocity scale ratio between the template pixel and the
                            nominal pixel width (default: 4)
      [1;32m-s[0m, [1;36m--step[0m [1;33mSTEP[0m       linear or log-linear spectral pixel width (default:
                            0.0001)
      [1;32m-l[0m, [1;36m--linear[0m          pixel sampling should be linear (default is log)
                            (default: False)
    