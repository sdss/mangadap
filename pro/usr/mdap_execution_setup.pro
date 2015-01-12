;===============================================================================
; STELLAR TEMPLATE LIBRARY FILES
;===============================================================================
;       The stellar template library files should be a list of 1D fits files,
;       and be associated with one of the following library keys (TODO: More to
;       be added):
;
;                                Spectral
;                 KEY    resolution (ang)       
;       -------------    ----------------
;           M11-MARCS                2.73
;          M11-STELIB                3.10
;
;       The fits files must have CRVAL1, CRPIX1, and CDELT1 keywords used to
;       define the wavelength coordinates of each pixel:
;
;               pix = (findgen(npixels+1))[1:*]
;               wave = (pix - CRPIX1) * CDELT1 + CRVAL1
;
;       The reference frame of the template wavelengths must also be defined as
;       either vacuum or air via the tpl_vacuum_wave vector.  Set the value to
;       be 1(true) if the wavelengths are in vacuum, 0(false) otherwise.  It is
;       expected that the DRP spectra are in vacuum wavelengths.  The DAP will
;       therefore use the IDL routine AIRTOVAC to convert the template
;       wavelengths to vacuum if tpl_vacuum_wave = 0.
;
;===============================================================================
;
;===============================================================================
; EMISSION-LINE PARAMETER FILES
;===============================================================================
;               File name with information regarding the emission lines to fit.
;               The format must be compatible with the module that performs the
;               fit of the emission lines.  If Gandalf or S-Gandalf are used,
;               the input file must be an ascii file with 9 columns as in the
;               following example (comments starts with ''\#''):
;
;       #  ID     CODE    WAVELENGTH   ACTION    LINE/  INTENSITY   V_g/i   sig_g/i   MODE
;       #                 [angstrom]  i/f/m/s  doublet                                f/tN
;       #   0     HeII       3203.15        m        l      1.000       0       10     t25
;       #   1     NeV        3345.81        m        l      1.000       0       10     t25
;           2     NeV        3425.81        m        l      1.000       0       10     t25
;           3     OII        3726.03        m        l      1.000       0       10     t25
;              
;               The columns are:
;
;               1. ID: Unique integer identifyer of the emission line.
;                       
;               2. CODE: Name of the element, read as a string.  Special
;               characters, such as '[', '.' are not permitted.
;                       
;               3. WAVELENGTH: Rest frame wavelength of the emission line to
;               fit.  WARNING: the name of the emission line in the DAP output
;               is: CODE+'_'+MDAP_STC(ROUND(wav), /integer)  
;                       
;               4. ACTION: Describes how the line should be treated.  Possible
;               values are:
;
;                       'i': ignore the line, as if the line were commented out.
;
;                       'f': fit the line and mask the line when fitting the
;                       stellar continuum.
;
;                       'm': mask the line when fitting the stellar continuum
;                       but do NOT fit the line itself
;
;                       's': defines a sky line that should be masked.  When
;                       masked, the wavelength of the line is NOT adjusted for
;                       the redshift of the object spectrum.
;                       
;               5. LINE:  Type of line, which can be either 'l' for a line or
;               'dN' for a doublet.  For doublets, the 'N' indicates the line ID
;               of the other line in the doublet.  The line to which the doublet
;               is tied should have the LINE='l'; for example, if emission line
;               with ID=4 has line='d3', then the emission line with ID=3 must
;               have LINE='l'.
;                       
;               6. INTENSITY:  Relative intensity of the gas emission (positive)
;               or absorption (negative) lines with respect to the doublet.
;               Therfore, this should most often be unity if LINE='l' and
;               indicate the ratio of line INTENSITY if LINE='dN'.
;                       
;               7. V_g/i: Guess for the velocity offset with respect the galaxy
;               systemic velocity.  TODO: This value is currently ignored by the
;               DAP!
;                       
;               8. sig_g/i. Guess for the velocity dispersion.  TODO: This value
;               is currently ignored by the DAP!
;                       
;               9. MODE.  Fitting mode for the line, which can be either 'f' to
;               fit the line independently or 'tN' to set both the velocity and
;               dispersion to be tied to a fitted line (MODE='f') with ID=N.
;               One can also force only the velocities to be tied using 'vN' or
;               only the velocity dispersions to be tied using 'sN'.
;
;       The reference frame of the emission-line wavelengths must also be
;       defined as either vacuum or air via the ems_vacuum_wave vector.  Set the
;       value to be 1(true) if the wavelengths are in vacuum, 0(false)
;       otherwise.  It is expected that the DRP spectra are in vacuum
;       wavelengths.  The DAP will therefore use the IDL routine AIRTOVAC to
;       convert the emission-line wavelengths to vacuum if ems_vacuum_wave = 0.
;
;===============================================================================
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; SPECTRAL INDEX PARAMETER FILES
;===============================================================================
;       Files must be ascii text with the following 9 columns, e.g.:
;
;       # ID Name        pass_0    pass_1   blcnt_0   blcnt_1   rdcnt_0   rdcnt_1 units
;       #                   ang       ang       ang       ang       ang       ang
;          0 D4000      0000.00   0000.00   3750.00   3950.00   4050.00   4250.00   ang 
;          1 CaII0p39   3899.50   4003.50   3806.50   3833.80   4020.70   4052.40   ang
;          2 HDeltaA    4083.50   4122.25   4041.60   4079.75   4128.50   4161.00   ang
;          3 HDeltaF    4091.00   4112.25   4057.25   4088.50   4114.75   4137.25   ang
;
;               1. Integer. Unique ID number of the absorption line feature.
;
;               2. String. Unique name of the absorption line feature. This will
;               define the name of the field in sctructure of the DAP results
;               (i.e. the name must begin with a letter, special characters like
;               comas or dots are not allowed).
;
;               3-4. Float (units: ang) Lower and upper value of the index
;               passband.  If these two number are both less than one (see
;               MDAP_INDEX_IS_BANDHEAD), the index is treated as a bandhead or
;               spectral break.
;
;               5-6. Float (units: ang) Lower and upper value of the index blue
;               pseudo-continuum.
;
;               7-8. Float (units: ang) Lower and upper value of the index red
;               pseudo-continuum.
;
;               9. String (accepted values are: ang or mag). Specifies the units
;               (ang or magnitudes) of the output.
;
;       The reference frame of the absorption-line wavelength parameters must
;       also be defined as either vacuum or air via the abs_vacuum_wave vector.
;       Set the value to be 1(true) if the wavelengths are in vacuum, 0(false)
;       otherwise.  It is expected that the DRP spectra are in vacuum
;       wavelengths.  The DAP will therefore use the IDL routine AIRTOVAC to
;       convert the absorption-line wavelength parameters to vacuum if
;       ems_vacuum_wave = 0.
;
;       NOTE: Indices will be measured only if their blue and red
;       pesudo-continua bandpasses are included in the considered wavelength
;       range. If not, their values are set to NaN, and their errors to 99 in
;       the final output file.
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; SPATIAL BINNING
;===============================================================================
;
;       There are currently four types of spatial binning allowed by the
;       DAP.  This is a list of the options and the parameters in the
;       BinPar structure that MUST accompany any usage of each binning
;       type.  In practice, if bin_par.type is not 'ALL', 'STON', or
;       'RADIAL' it is assumed to be 'NONE'.
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
;           bin_par.type = 'STON'
;               Spectra are binned, using the Voronoi binning scheme, to
;               a minimum S/N level. In addition to setting
;               bin_par.v_register and bin_par.optimal_weighting (see
;               above), this bin type also requires:
;
;                   bin_par.ston
;                       Minimum S/N level.
;
;           bin_par.type = 'RADIAL'
;               Spectra are binned radially according to a provided
;               planar projection.  In addition to setting
;               bin_par.v_register and bin_par.optimal_weighting (see
;               above under bin_par.type='ALL'), this bin type also
;               requires:
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
;===============================================================================
;===============================================================================
;
;===============================================================================
; ANALYSIS STEPS
;===============================================================================
;       The currently available analysis steps are:
;
;           'stellar-cont'
;               Fit the stellar continuum using pPXF.  Requires the
;               template library.  Will mask emission lines, if a set of
;               GANDALF-style emission-line parameters are provided.
;               Determines the optimal template mix and the stellar
;               kinematics.
;
;           'star+gas'
;               This step uses GANDALF to fit the spectrum, which will
;               re-optimize the template mix, but does not adjust
;               the stellar kinematics.  Determines the emission-line
;               kinematics, intensities, fluxes, and equivalent widths.
;
;           'emission-line'
;               Fits an emission-line-only spectrum assuming a zero
;               baseline.  If either 'stellar-cont' or 'star+gas' have
;               been applied, the best-fitting stellar continuum (where
;               the GANDALF results take preference) will be subtracted
;               from the input spectra before fitting the
;               emission-line-only model.  If this has not been
;               performed, **no continuum is subtracted from the
;               spectra!**  This fit, therefore, does not require a
;               template library to be provided.  Currently, the lines
;               fit are hard-coded to be:
;
;               #   Line    Rest Wave (air)
;                   OII     3726.03
;                   Hb      4861.32
;                   OIII    4958.83
;                   OIII    5006.77
;                   NII     6547.96
;                   Ha      6562.80
;                   NII     6583.34
;                   SII     6716.31
;                   SII     6730.68
;
;               The fit is done twice using two different contributed
;               codes, one from Enci Wang and another from Francesco
;               Belfiore.  For these lines only, this step determines
;               the emission-line kinematics, intensities, fluxes, and
;               equivalent widths.
;
;           'abs-indices'
;               Calculates the absorption-line indices.  Requires an
;               absorption-line parameter set and a fit to the stellar
;               continuum.  If an emission-line model has been fit (with
;               the one from GANDALF taking precedence), it will be
;               subtracted from the galaxy spectra before performing the
;               measurements.  Other steps applied before performing the
;               index measurements is to replace aberrant pixels and to
;               match the spectral resolution to that of the index
;               system.  Index measurements performed on the
;               best-fitting template and its broadened version are used
;               to correct the measurements made for the data in a
;               differential way.  This step provides index measurements
;               and their errors.
;
; NOTE: The current list of fits must occur sequentially, except the
; emission-line only fit, which can be done without any continuum fit,
; although this is not recommended!  Based on the requested analysis
; steps, the code has logic in place that will determine what other
; steps need to be performed, and will do so after warning the user.
;===============================================================================
;===============================================================================
;
;===============================================================================
; ANALYSIS PARAMETERS
;===============================================================================
;
;       Some of the analysis steps have parameters that can be set to
;       specify certain aspects of their execution.  These parameters
;       are collected into a single structure called AnalysisPar and
;       defined in mdap_define_analysis_par.pro.  The currently
;       available parameters are:
;
;           AnalysisPar.moments
;               Number of *stellar* kinematics to fit.  The number of
;               emission-line kinematic moments is currently fixed at 2
;               (v, sigma).
;
;           AnalysisPar.degree
;               Degree of the *additive* polynomial used by pPXF when
;               fitting the stellar continuum.
;
;           AnalysisPar.mdegree
;               Degree of the *multiplicative* polynomial used by pPXF
;               and GANDALF when fitting the stellar continuum and full
;               spectrum, respectively.
;
;           AnalysisPar.reddening_order
;               The order of the reddening fit by GANDALF:  0 means no
;               reddening; 1 means fit a single dust-screen model; 2
;               means fit a dust-screen model plus a nebular only
;               component.  NOTE, if the reddening is fit, the value of
;               MDEGREE is ignored because fitting them both is
;               degenerate.
;
;           AnalysisPar.reddening
;               A two-element array used to set the initial guess for
;               the reddening coefficients.  The provided values must
;               match the input reddening_order, but the array must
;               always contain two elements.
;
;           AnalysisPar.zero_instr_disp
;               A flag to force GANDALF to ignore the instrumental
;               dispersion (set it to 0).  Flag is 0-false,1-true;
;               default is false.
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; ANALYSIS PRIORS
;===============================================================================
;
;       When executing a certain execution plan, you may want to use
;       results from a previous run of the DAP as priors for another
;       execution plan.  You can do this through the definition of the
;       analysis_prior.  This is a string array that contains either the
;       file name of the DAP fits file to use (must use the full path or
;       the path from the execution directory) or the index number (in
;       string format) of the execution plan in the current list to use
;       as the prior.
;        
;       For example, for two plans where the first has no prior and the
;       second uses the first as a prior, you would define:
;
;           analysis_prior[0] = ''
;           analysis_prior[1] = '0'
;
;===============================================================================
;===============================================================================
;
; Setup some necessary execution variables for the MANGA DAP

; Signifier is what will be reported to the log file as the configuration file
; for a run of MaNGA_DAP

; dapsrc is an optional input to define the DAP source path instead of
; using environmental varaibles.

PRO MDAP_EXECUTION_SETUP, $
        tpl_library_keys, template_libraries, tpl_vacuum_wave, ems_line_keys, $
        emission_line_parameters, ems_vacuum_wave, abs_line_keys, absorption_line_parameters, $
        abs_vacuum_wave, signifier, bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, $
        threshold_ston_analysis, analysis, tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
        analysis_par, analysis_prior, overwrite_flag, dapsrc=dapsrc, $
        save_intermediate_steps=save_intermediate_steps, $
        remove_null_templates=remove_null_templates, external_library=external_library


        ; Define the DAP source path
        if n_elements(dapsrc) eq 0 then $
            dapsrc = getenv('MANGADAP_DIR')

        ;-----------------------------------------------------------------------
        ; Location for the output files
        ;
        ;       For a file name of 
        ;               manga-7443-12701-LOGCUBE.fits
        ;
        ;       Output files will be placed in
        ;               output_root_dir+'/'+manga-7443-12701-LOGCUBE
        ;
        ; output_root_dir=getenv('MANGA_SPECTRO_ANALYSIS')+'/'+getenv('MANGADAP_VER')
        ;
        ;   OBSOLETE! Output path now either defined by the default or
        ;   dappath set when calling MANGA_DAP

        ;-----------------------------------------------------------------------
        ; List of galaxies to process and their properties ASCII file, one
        ; line per galaxy, format:
        ;    Col 1. [string] Name (and location) of the fits file, WITHOUT the
        ;                    .fits extension
        ;    Col 2. [float]  Galaxy mean velocity (km/sec)
        ;    Col 3. [float]  Galaxy mean velocity dispersion (km/sec)
        ;    Col 4. [float]  Galaxy mean ellipticity
        ;    Col 5. [float]  Galaxy mean Position Angle
        ;    Col 6. [float]  MANGA Bundle size
        ;    Col 7. [float]  Galaxy effective radius
        ;    Col 8. [string] Dataformat. Valid entries are: CUBE or RSS
        ;
        ; total_filelist=output_root_dir+'/manga_dap_table.inp' 
        ;
        ;   OBSOLETE! File now not needed if using a par file or defined
        ;   using the inptbl when calling MANGA_DAP

        ;-----------------------------------------------------------------------
        ; Flag to save intermediate steps.
        ; TODO: This is no longer used!

        ; save_intermediate_steps = 0 

        ;-----------------------------------------------------------------------
        ; Remove templates with zero weights in one fit from use in another fit.
        ; TODO: Currently not implemented.  Will include this as an option when
        ; applying priors.

        ; remove_null_templates = 1

        ;-----------------------------------------------------------------------
        ; Path to a library of fortran or C codes to be used.  If commented,
        ; internal IDL procedures are used.

        ; external_library=getenv('MANGADAP_DIR')+'/external/F90_32/'
        ; external_library=getenv('MANGADAP_DIR')+'/external/F90_64/'
        external_library=dapsrc+'/external/F90_64/'

        ;-----------------------------------------------------------------------
        ; Define the set of template libraries.  The format expected for these
        ; files is described above.
        ntpl_libraries = 2
        tpl_library_keys = strarr(ntpl_libraries)
        template_libraries = strarr(ntpl_libraries)
        tpl_vacuum_wave = intarr(ntpl_libraries)

        tpl_library_keys[0] = 'M11-MARCS'
        template_libraries[0] = dapsrc+'/external/templates_m11_marcs/*_s.fits'
        tpl_vacuum_wave[0] = 0

        tpl_library_keys[1] = 'M11-STELIB'
        template_libraries[1] = dapsrc+'/external/templates_m11_stelib/*_s.fits'
        tpl_vacuum_wave[1] = 0

        ;-----------------------------------------------------------------------
        ; Define the set of emission-line parameter files.  The format expected
        ; for these files is described above.
        neml_files = 1
;       neml_files = 3
        ems_line_keys = strarr(neml_files)
        emission_line_parameters = strarr(neml_files)
        ems_vacuum_wave = intarr(neml_files)

        ems_line_keys[0] = 'STANDARD'
        emission_line_parameters[0] = dapsrc+'/external/manga_emission_line_list_nominal.par'
        ems_vacuum_wave[0] = 1

;       ems_line_keys[1] = 'NODOUBLETS'
;       emission_line_parameters[1] = dapsrc + $
;                       '/external/emission_lines_setup_with_Balmer_decrement_no_doublets'
;       ems_vacuum_wave[1] = 0

;       ems_line_keys[2] = 'RESIDUAL'
;       emission_line_parameters[2] = dapsrc + $
;                       '/external/emission_lines_setup_with_Balmer_decrement_residuals'
;       ems_vacuum_wave[2] = 0

        ;-----------------------------------------------------------------------
        ; Define the set of absorption-line parameter files.  The format expected
        ; for these files is described above.
        nabs_files = 1
        abs_line_keys = strarr(nabs_files)
        absorption_line_parameters = strarr(nabs_files)
        abs_vacuum_wave = intarr(nabs_files)

        abs_line_keys[0] = 'LICK'
        absorption_line_parameters[0] = dapsrc+'/external/absorption_line_indices_definition.dat'
        abs_vacuum_wave[0] = 0

        ;=======================================================================
        ; DEFINITION OF EXECUTION PROCEDURES

        ; Define a string used to signify this file in the header of the DAP output file(s)
;       cd, current=directory
;       signifier = directory+'/mdap_setup.pro'
        signifier = dapsrc+'/pro/usr/mdap_execution_setup.pro'

        ;-----------------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.

        niter = 2                                       ; Number of ExecutionPlans to produce
;       niter = 4                                       ; Number of ExecutionPlans to produce

        bin_par_def = MDAP_DEFINE_BIN_PAR()             ; Define the BinPar structure
        bin_par = replicate( bin_par_def, niter)        ; Create the array of BinPar structures

        w_range_sn = dblarr(niter, 2)                   ; Wavelength range for S/N calculation
        threshold_ston_bin = dblarr(niter)              ; Threshold S/N to include spectrum in bin

        w_range_analysis = dblarr(niter, 2)             ; Wavelength range for the analysis
        threshold_ston_analysis = dblarr(niter)         ; Threshold S/N to analyze spectrum

        max_analysis_blocks = 4                         ; Maximum number of analysis blocks
        analysis = strarr(niter, max_analysis_blocks)   ; Analysis steps to apply

        tpl_lib_analysis = intarr(niter)                ; INDEX of template library to use
        ems_par_analysis = intarr(niter)                ; INDEX of emission-line parameter file
        abs_par_analysis = intarr(niter)                ; INDEX of absorption-line parameter file

        analysis_par_def = MDAP_DEFINE_ANALYSIS_PAR()   ; Define the AnalysisPar structure
        analysis_par = replicate( analysis_par_def, niter)  ; Create array of AnalysisPar structs

        analysis_prior = strarr(niter)                  ; Prior information used for analysis

        overwrite_flag = intarr(niter)                  ; Flag to overwrite any existing output file

        ;-----------------------------------------------------------------------
        ; For each iteration:
        
        ; Define the necessary elements of the BinPar structure.  The
        ; available binning types and their required parameters are
        ; provided above.


;-------------------------------------------------------------------------------
; To fit things similarly to the PDAP:
;       bin_par[0:2].type = 'STON'
;       bin_par[0].ston = 40.0d
;       bin_par[1].ston = 25.0d
;       bin_par[2].ston = 15.0d

;       bin_par[3].type = 'RADIAL'
;       bin_par[3].v_register = 1
;       bin_par[3].nr = 6
;       bin_par[3].rs = 1.0
;       bin_par[3].rlog = 1

;       bin_par[*].optimal_weighting = 1

;       w_range_sn[0,*] = [5560.00, 6942.00]
;       w_range_sn[1,*] = [5560.00, 6942.00]
;       w_range_sn[2,*] = [5560.00, 6942.00]
;       w_range_sn[3,*] = [5560.00, 6942.00]

;       threshold_ston_bin[*] = -300.0d

;       w_range_analysis[0,*] = [3650.,10300.]
;       w_range_analysis[1,*] = [3650.,7500.] 
;       w_range_analysis[2,*] = [3600.,10220.]
;       w_range_analysis[3,*] = [3700.,9980.]

;       threshold_ston_analysis[*] = 0.0d

;       analysis[*,0] = 'stellar-cont'
;       analysis[*,1] = 'star+gas'
;       analysis[*,2] = 'emission-line'
;       analysis[*,3] = 'abs-indices'

;       tpl_lib_analysis[*] = 0
;       ems_par_analysis[*] = 0
;       abs_par_analysis[*] = 0

;       analysis_par[0:2].moments = 4
;       analysis_par[3].moments = 2
;       analysis_par[0].degree = -1
;       analysis_par[1:2].degree = 3
;       analysis_Par[3].degree = -1
;       analysis_par[*].mdegree = 6
;       analysis_par[0:1].reddening_order = 0
;       analysis_par[3].reddening_order = 0
;       analysis_par[2].reddening_order = 2
;       analysis_par[2].reddening[*] = [0.01,0.01]

;       analysis_prior[0] = ''
;       analysis_prior[1] = '0'
;       analysis_prior[2] = '1'
;       analysis_prior[3] = '2'

;       overwrite_flag[*] = 0

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
; To fit just the high S/N case (with both MARCS and STELIB):
        bin_par[*].type = 'STON'
        bin_par[*].optimal_weighting = 1        ; Otherwise uniform weighting
        bin_par[*].ston = 40.0d

        w_range_sn[0,*] = [5560.00, 6942.00]
        w_range_sn[1,*] = [5560.00, 6942.00]
        threshold_ston_bin[*] = -300.0d

        w_range_analysis[0,*] = [3650.,10300.] 
        w_range_analysis[1,*] = [3650.,10300.] 
        threshold_ston_analysis[*] = 0.0d

        analysis[*,0] = 'stellar-cont'
        analysis[*,1] = 'star+gas'
        analysis[*,2] = 'emission-line'
        analysis[*,3] = 'abs-indices'

        tpl_lib_analysis[0] = 0
        tpl_lib_analysis[1] = 1
        ems_par_analysis[*] = 0
        abs_par_analysis[*] = 0

        analysis_par[*].moments = 4
        analysis_par[*].degree = -1
        analysis_par[*].mdegree = 6
        analysis_par[*].reddening_order = 0
        analysis_par[*].zero_instr_disp = 1

        analysis_prior[*] = ''      ; No priors

        overwrite_flag[*] = 1

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


; Another set of posibilites, with more comments as to what you're setting
;       bin_par[0].type = 'STON'
;       bin_par[0].optimal_weighting = 1        ; Otherwise uniform weighting
;       bin_par[0].ston = 40.0d
        ;   leave everything else as default (no velocity registration)

        ; Try RADIAL using the results from the first ExecutionPlan to
        ; velocity register the data -> set v_register to true here and
        ; add the prior below.
;       bin_par[1].type = 'RADIAL'
;       bin_par[1].v_register = 1
;       bin_par[1].optimal_weighting = 1
;       bin_par[1].nr = 10
;       bin_par[1].rlog = 1
        ;   leave everything else as default

        ; Define the wavelength range over which to calculate the mean S/N per pixel
;       w_range_sn[0,*] = [5560.00, 6942.00]
;       w_range_sn[1,*] = [5560.00, 6942.00]

        ; Define the S/N threshold to include spectrum in any bin
;       threshold_ston_bin[*] = -300.0d

        ; Define the wavelength range over which to perform ALL analyses
;       w_range_analysis[0,*] = [3650.,10300.] 
;       w_range_analysis[1,*] = [3650.,10300.] 

        ; Define the S/N threshold to perform analysis
;       threshold_ston_analysis[*] = 0.0d

        ; Set the list of analyses to perform.  The available analysis steps are
        ; listed above.

;       analysis[*,0] = 'stellar-cont'
;       analysis[*,1] = 'star+gas'
;       analysis[*,2] = 'emission-line'
;       analysis[*,3] = 'abs-indices'

        ; Set the index of the template library to use for the analysis
        ; TODO: Change this to use the library key?
;       tpl_lib_analysis[*] = 0

        ; Set the index of the emission-line parameter set to use
;       ems_par_analysis[*] = 0

        ; Set the index of the absorption-line parameter set to use
;       abs_par_analysis[*] = 0

        ; Set additional parameters needed by the analysis modules
        ; The reddening order can be 0, 1, or 2
        ; TODO: Allow for oversample?
        ; IF NOT SET HERE, the default values are:
        ;       moments=2, degree=-1, mdegree=-1, reddening_order=0
;       analysis_par[*].moments = 4
;       analysis_par[*].degree = -1
;       analysis_par[*].mdegree = 6
;       analysis_par[*].reddening_order = 0
        ; analysis_par[0].reddening[*] = [0.01,0.01]

        ; Analysis priors, see description above.
;       analysis_prior[0] = ''      ; No prior for the first plan
;       analysis_prior[1] = '0'     ; Use the results from the first plan as a prior on the second

        ; Set a flag to overwrite existing output: 1-yes, 0-no
;       overwrite_flag[0] = 1
;       overwrite_flag[0] = 1
;       overwrite_flag[1] = 1

        ;=======================================================================
        ;=======================================================================

END

