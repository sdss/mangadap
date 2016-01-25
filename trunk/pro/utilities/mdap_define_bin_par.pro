;+
; NAME:
;       MDAP_DEFINE_BIN_PAR
;
; PURPOSE:
;       Define the binning parameter structure and its default values
;       for use.  This structure holds all binning parameters including
;       the type of binning to perform.
;
;       There are currently four types of spatial binning allowed by the
;       DAP.  This is a list of the options and the parameters in the
;       BinPar structure that MUST accompany any usage of each binning
;       type.  In practice, bin_par.type is assumed to be 'NONE' if it
;       is not 'ALL', 'STON', or 'RADIAL'.
;
;           bin_par.type='NONE'
;               No binning is performed, all spectra are analyzed
;               (assuming they meet the analysis criteria).  No
;               parameters are required for this binning type.
;
;           bin_par.type='ALL'
;               All spectra are combined into a single bin for analysis.
;               The parameters required for this binning type are:
;
;                   bin_par.v_register
;                       Use a prior fit of the kinematics to velocity
;                       register the spectra to the median redshift
;                       before binning
;
;                   TODO: Allow for a bin_par.v_tweak keyword that will
;                   use a function to tweak the velocity registration.
;
;                   bin_par.optimal_weighting
;                       Weight spectra by S/N^2 during their combination
;                       (1-yes,0-no)
;
;                   bin_par.noise_calib
;                       Apply a calibration of the noise vector in the
;                       binned spectrum.  Options are currently:
;                           
;                           0 - (default) Do not apply
;                           1 - N_calib = N_nominal * (1 + 1.6*log10(N_bin))
;
;                       where N_bin is the number of binned spectra.
;                       NOTE: Option 1 cannot be combined with optimal
;                       weighting!
;
;           bin_par.type = 'STON'
;               Spectra are binned, using the Voronoi binning scheme, to
;               a minimum S/N level. In addition to setting
;               bin_par.v_register, bin_par.optimal_weighting, and
;               bin_par.noise_calib (see above), this bin type also
;               requires:
;
;                   bin_par.ston
;                       Minimum S/N level.
;
;           bin_par.type = 'RADIAL'
;               Spectra are binned radially according to a provided
;               planar projection.  In addition to setting
;               bin_par.v_register, bin_par.optimal_weighting,
;               bin_par.noise_calib (see above under
;               bin_par.type='ALL'), this bin type also requires:
;
;                   bin_par.cx, bin_par.cy:
;                       On-sky X,Y position to use for R=0.  (0,0 by
;                       default)
;
;                   bin_par.pa:
;                       Position angle from North (0 deg) through east
;                       (90 deg) of the MAJOR axis of a constant radius
;                       ellipse. (0 by default)
;
;                   bin_par.ell:
;                       Ellipticity (1-b/a) of a constant radius
;                       ellipse. (0 by default)
;
;                   bin_par.rs:
;                       Starting radius of the first bin (0 by default)
;
;                   bin_par.re:
;                       Ending radius of the last bin (largest R by
;                       default)
;
;                   bin_par.nr:
;                       Number of radial bins (5 by default)
;
;                   bin_par.rlog:
;                       Flag to logarithmically step the radial bins
;                       (1-true;0-false -- 0 by default)
;
;                   bin_par.rscale:
;                       Value by which to scale the radius in arcsec
;                       (1.0 by default)
;
;               The values of bin_par.pa and bin_par.ell WILL BE SET
;               AUTOMATICALLY using the ell and pa parameters read by
;               MANGA_DAP.  bin_par.rscale can be set to Reff
;               automatically by setting bin_par.rscale = -1.0.
;               
;               TODO: Add more options for rscale...
;
;       TODO: Other binning options maybe added later...
;
;       TODO: Include a 'redshift_prior' parameter that points to a file
;       or ExecutionPlan iteration to use to de-redshift the spectra
;       before binning them.  For 'NONE', the spectra are just
;       deredshifted.
;
; CALLING SEQUENCE:
;       result = MDAP_DEFINE_BIN_PAR()
;
; REVISION HISTORY:
;       04 Dec 2014: Original implementation by K. Westfall.  bin_type
;                    and bin_weight_by_sn2 are no longer parameters in
;                    the ExecutionPlan structure.  Instead they're held
;                    here.
;       16 Mar 2015: (KBW) Added noise_calib keyword for the S/N
;                          calibration of binned spectra that accounts
;                          for covariance.
;-
;-----------------------------------------------------------------------

FUNCTION MDAP_DEFINE_BIN_PAR
        bin_par = { BinPar, type:'', v_register:0, optimal_weighting:0, noise_calib:0, ston:0.0d, $
                            cx:0.0d, cy:0.0d, pa:0.0d, ell:0.0d, rs:0.0d, re:-1.0d, nr:5, rlog:0, $
                            rscale:1.0d }

        return, bin_par
END

