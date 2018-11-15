drpf = DRPFits(obs['plate'], obs['ifudesign'], obs['mode'], drpver=_drpver,
                 redux_path=redux_path, directory_path=directory_path, read=True)

rdxqa = ReductionAssessment('SNRR', drpf, pa=obs['pa'], ell=obs['ell'], clobber=False)

binned_spectra = SpatiallyBinnedSpectra('VOR30', drpf, rdxqa, reff=obs['reff'], clobber=False)

stellar_continuum = StellarContinuumModel('GAU-MILESSP', binned_spectra, guess_vel=obs['vel'], guess_sig=obs['vdisp'],
                                          clobber=False)
