;+
; NAME:
;       MDAP_BUILD_EXECUTION_PLANS
;
; PURPOSE:
;       Populate a set of ExecutionPlan structures defining how the DAP is
;       executed for a given DRP-produced fits file.
;
; CALLING SEQUENCE:
;       MDAP_BUILD_EXECUTION_PLANS, n_tpl, n_ems, n_abs, bin_par, ell, pa, Reff, w_range_sn, $
;                                   threshold_ston_bin, w_range_analysis, threshold_ston_analysis, $
;                                   analysis, analysis_par, analysis_prior, tpl_lib_analysis, $
;                                   ems_par_analysis, abs_par_analysis, overwrite_flag, file_root, $
;                                   execution_plan, version=version
;
; INPUTS:
;       n_tpl integer
;               Number of available template libraries.  This number is checked
;               against the selected template library for the analysis
;               (tpl_lib_analysis) for each plan, if a library is required.
;
;       n_ems integer
;               Number of available emission-line parameter sets.  This number
;               is checked against the selected parameter set for the analysis
;               (ems_par_analysis) for each plan, if a parameter set is
;               required.
;
;       n_abs integer
;               Number of available absorption-line parameter sets.  This number
;               is checked against the selected parameter set for the analysis
;               (abs_par_analysis) for each plan, if a parameter set is
;               required.
;
;       bin_par BinPar[P]
;               Structure containing the parameters used to perform the
;               spatial binning.  See MDAP_DEFINE_BIN_PAR() for a
;               description of the structure.
;
;       ell double
;               Ellipticity (1-b/a) of the galaxy.  Copied to
;               bin_par.ell for any of the P BinPar structures with type
;               == 'RADIAL'.
;
;       pa double
;               Position angle of the galaxy.  Copied to bin_par.pa for
;               any of the P BinPar structures with type == 'RADIAL'.
;
;       Reff double
;               Effective radius of the galaxy.  Copied to
;               bin_par.rscale for any of the P BinPar structures with
;               type == 'RADIAL' and rscale le 0.0.
;                       
;       w_range_sn dblarr[P][2]
;               Wavelength range used to calculate the S/N of each spectrum,
;               performed using each individual spectrum.  This is always
;               calculated, regardless of the binning type.
;
;       threshold_ston_bin dblarr[P]
;               The S/N threshold for the inclusion of any DRP spectrum in the
;               binning process.
;
;       w_range_analysis dblarr[P][2]
;               The wavelength range to use during the analysis.
;
;       threshold_ston_analysis dblarr[P]
;               The S/N threshold for the inclusion of a spectrum in the
;               analysis.
;
;       analysis strarr[P][]
;               List of analyses to perform.  Valid values are:
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
;       analysis_par AnalaysisPar[P]
;               An array of structures that keeps the parameters used by
;               the analysis steps.  See MDAP_DEFINE_ANALYSIS_PAR().
;
;       analysis_prior strarr[P]
;               A string list of DAP files to use as a prior during each
;               execution plan.  On input, these are expected to contain
;               either the name of the DAP file directly (with either
;               the full path or the path within the calling directory)
;               or a string representation of an index of the plan to
;               use for the prior.  (**THIS FUNCTION**) replaces the
;               index numbered priors with the full file name, which is
;               what is included in the ExecutionPlan structure.
;
;       tpl_lib_analysis intarr[P]
;               Index of template library to use for each of the P execution
;               plans.  -1 if no template library is to be used during the
;               analysis.
;
;       ems_par_analysis intarr[P]
;               Index of emission-line parameter set to use for each of the P
;               execution plans.  -1 if no emission-line parameter set is to be
;               used during the analysis.
;
;       abs_par_analysis intarr[P]
;               Index of absorption-line parameter set to use for each of the P
;               execution plans.  -1 if no absorption-line parameter set is to
;               be used during the analysis.
;
;       overwrite_flag intarr[P]
;               Flag to overwrite any existing output file with the exact same
;               execution plan.  TODO: Should this actually just be a flag that
;               forces analyses to be redone?
;
;       file_root string
;               Root name for the output file of each plan.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       execution_plan ExecutionPlan[P]
;               Array of P structures containing the execution plans.
;
; OPTIONAL OUTPUT:
;       version string
;               Module version. If requested, the module is not executed and only
;               the version flag is returned
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       10 Oct 2014: Original implementation by K. Westfall (KBW)
;       13 Oct 2014: (KBW) Include bin_weight_by_sn2
;       15 Oct 2014: (KBW) Added analysis_extra
;       07 Nov 2014: (KBW) Checks that analysis steps are performed in
;                          sequential order is now done using
;                          MDAP_ANALYSIS_BLOCKS_TO_PERFORM(), which takes into
;                          account the blocks that have already been completed.
;       04 Dec 2014: (KBW) Changes to allow for RADIAL binning type,
;                          which is mostly a change to bin_par.
;       05 Dec 2014: (KBW) Include analysis prior in execution plan,
;                          moved MDAP_GENERATE_OUTPUT_FILE_NAMES from
;                          MANGA_DAP to here and converted it to a
;                          function, added file_root to input, added
;                          versioning.
;       11 Dec 2014: (KBW) Allow for the emission-line only fitting.
;       16 Mar 2015: (KBW) Allow noise_calib in bin_par.  Force
;                          optimal_weighting = 0 if applying the
;                          calibration.
;       22 Mar 2015: (KBW) Include a check of the number of pPXF
;                          moments.  Automatically set the number of
;                          moments to 2 if input is set to 0.
;-
;------------------------------------------------------------------------------

; TODO: Eventually this will generate the files names AND check them against
; existing files.  Any execution plan that matches IN DETAIL the execution plan
; of an existing file (apart from the analyses performed) will be used to input
; already completed analyses, unless the overwrite flag is flipped.   This check
; will *include* checks of the module versions used to generate the existing
; analysis.
;
; Given the stuff in mdap_analysis_block_logic.pro, I'm not sure the
; above is true anymore...

;-----------------------------------------------------------------------
; Below used to be MDAP_GENERATE_OUTPUT_FILE_NAMES
;
FUNCTION MDAP_OUTPUT_FILE_NAME, $
                file_root, execution_plan, i
        return, (file_root + 'BIN-' + execution_plan.bin_par.type + '-' + $
                string(i+1,format='(I03)') + '.fits')
END
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
FUNCTION MDAP_SET_EXECUTION_PLAN_ANALYSIS, $
                analysis, n_analyses

        n_steps = n_elements(analysis)          ; Number of analysis steps

        analysis_ = intarr(n_analyses)          ; Flag array used to toggle the analysis steps
        if n_steps eq 0 then $
            return, analysis_                   ; All analysis steps set to zero

        if n_steps gt n_analyses then $
            message, 'Too many analysis steps provided!'

        ; Set the steps
        for i=0,n_steps-1 do begin
            if analysis[i] eq 'stellar-cont' then begin
                analysis_[0] = 1                ; Perform the stellar continuum fitting
                continue
            endif
            if analysis[i] eq 'star+gas' then begin
                analysis_[1] = 1                ; Perform the "simultaneous" star+gas fit
                continue
            endif
            if analysis[i] eq 'emission-line' then begin
                analysis_[2] = 1                ; Perform the emission-line-only fitting
                continue
            endif
            if analysis[i] eq 'abs-indices' then begin
                analysis_[3] = 1                ; Perform the absorption index measurements
                continue
            endif
            
            if strlen(analysis[i]) ne 0 then $
                message, 'Unknown analysis step: ', analysis[i]
        endfor

        ; TODO: So far all steps must occur serially.  True?  Can
        ;       abs-indices be done without doing star+gas?

        ; TODO: This is now checked in MDAP_ANALYSIS_BLOCKS_TO_PERFORM()

;       analysis_out = analysis_
;       for i=2,0,-1 do begin
;           if analysis_out[i] eq 1 then begin
;               for j=i-1,0,-1 do $
;                   analysis_out[j] = 1
;           endif
;       endfor
;
;       ; Notify the user that changes were made
;       indx = where(analysis_out - analysis_ ne 0)
;       if indx[0] ne -1 then begin
;           print, 'WARNING: Analysis steps are sequential.  Changed execution plan to ' + $
;                  'compensate.  Previous analyses will not be redone unless it was chosen to ' + $
;                  'overwrite the previous analysis.'
;       endif
;
;       return, analysis_out
        return, analysis_
END

PRO MDAP_CHECK_EXECUTION_PLAN, $
                n_tpl, n_ems, n_abs, execution_plan

        ; Check bin type
;       if execution_plan.bin_par.type ne 'NONE' and execution_plan.bin_par.type ne 'ALL' and $
;          execution_plan.bin_par.type ne 'STON' and execution_plan.bin_par.type ne 'RADIAL' $
        if execution_plan.bin_par.type ne 'NONE' && execution_plan.bin_par.type ne 'ALL' && $
           execution_plan.bin_par.type ne 'STON' && execution_plan.bin_par.type ne 'RADIAL' $
           then begin 
            message, 'Unknown binning type:' + execution_plan.bin_par.type
        endif

        ; If there is no prior, cannot velocity register the spectra
;       if execution_plan.bin_par.v_register eq 1 $
;          and strlen(execution_plan.analysis_prior) eq 0 then begin
        if execution_plan.bin_par.v_register eq 1 $
           && strlen(execution_plan.analysis_prior) eq 0 then begin
            message, 'Cannot velocity register spectra for binning without an analysis prior!'
        endif

;       ; Check bin parameters
;       if (execution_plan.bin_type eq 'STON' and n_elements(execution_plan.bin_par) ne 1) or $
;          (execution_plan.bin_type eq 'RAD' and n_elements(execution_plan.bin_par) ne 1) then begin
;           message, 'STON and RAD bin types must have one binning parameter (S/N and Bin width,'+ $
;                  ' respectively)'
;           
;       endif

        ; Check the noise_calib and optimal_weighting parameters are compatible
        if execution_plan.bin_par.optimal_weighting eq 1 $
           && execution_plan.bin_par.noise_calib ne 0 then begin
            print, 'WARNING: Cannot optimally weight spectra and apply the noise calibration.  '+$
                   'Ignoring optimal weighting.'
            execution_plan.bin_par.optimal_weighting = 0
        endif

        ; Analysis steps are check using MDAP_SET_EXECUTION_PLAN_ANALYSIS

        ; Same check as is done in pPXF
        if execution_plan.analysis[0] eq 1 $
           && total(execution_plan.analysis_par.moments eq [0,2,4,6]) eq 0 then $
            message, 'pPXF moments should be 0 (use default=2), 2, 4 or 6'

        ; Set the default number of moments
        if execution_plan.analysis[0] eq 1 && execution_plan.analysis_par.moments eq 0 then $
            execution_plan.analysis_par.moments = 2
       
        ; Check that the selected template library exists
        if execution_plan.tpl_lib ne -1 && execution_plan.tpl_lib ge n_tpl then $
            message, 'No template library index: ', execution_plan.tpl_lib
            
        ; Check that the selected emission-line parameter set exists
        if execution_plan.ems_par ne -1 && execution_plan.ems_par ge n_ems then $
            message, 'No emission-line parameter set: ', execution_plan.ems_par
            
        ; Check that the selected absorption-line parameter set exists
        if execution_plan.abs_par ne -1 && execution_plan.abs_par ge n_abs then $
            message, 'No absorption-line parameter set: ', execution_plan.abs_par


        ; If analysis parameters are not needed, then turn them off
        ;
        ; Absorption-line parameters only needed for that analysis
        if execution_plan.analysis[3] eq 0 then $
            execution_plan.abs_par = -1
        ; Emission-line parameters used for both stellar-continuum and star+gas analysis
        if execution_plan.analysis[1] eq 0 && execution_plan.analysis[0] eq 0 then $
            execution_plan.ems_par = -1
        ; Template spectra can be used for all analyses (what to do with [2]?)
        if execution_plan.analysis[3] eq 0 && execution_plan.analysis[2] eq 0 && $
           execution_plan.analysis[1] eq 0 && execution_plan.analysis[0] eq 0 then begin
            execution_plan.tpl_lib = -1
        endif

END

; TODO: A better "analysis structure"?
;           .name
;           .requested
;           .perform
;           .dependencies

PRO MDAP_BUILD_EXECUTION_PLANS, $
                n_tpl, n_ems, n_abs, bin_par, ell, pa, Reff, w_range_sn, threshold_ston_bin, $
                w_range_analysis, threshold_ston_analysis, analysis, analysis_par, analysis_prior, $
                tpl_lib_analysis, ems_par_analysis, abs_par_analysis, overwrite_flag, file_root, $
                execution_plan, version=version

        version_module = '0.3'                          ; Version number
        if n_elements(version) ne 0 then begin          ; set version and return
            version = version_module
            return
        endif

        n_plans = n_elements(bin_par)                   ; Number of plans to execute

        bin_par_def = MDAP_DEFINE_BIN_PAR()
        analysis_par_def = MDAP_DEFINE_ANALYSIS_PAR()

        n_analyses = 4                                  ; Number of analyses that can be selected

        ; Instantiate the ExecutionPlan structure
        ; TODO: Make analysis_extra a pointer?
        execution_plan = replicate( { ExecutionPlan, bin_par:bin_par_def, wave_range_sn:dblarr(2), $
                                                     threshold_ston_bin:0.0d, $
                                                     wave_range_analysis:dblarr(2), $
                                                     threshold_ston_analysis:0.0d, $
                                                     analysis:intarr(n_analyses), $
                                                     analysis_par:analysis_par_def, $
                                                     analysis_prior:'', tpl_lib:0, ems_par:0, $
                                                     abs_par:0, overwrite:0, ofile:''}, n_plans)

        for i=0,n_plans-1 do begin

            execution_plan[i].bin_par = bin_par[i]                      ; Binning parameters
            if bin_par[i].type eq 'RADIAL' then begin
                execution_plan[i].bin_par.ell = ell
                execution_plan[i].bin_par.pa = pa
                if execution_plan[i].bin_par.rscale le 0.0 then $
                    execution_plan[i].bin_par.rscale = Reff
                ; TODO: Check if prior exists.  Create warning if it does not?
            endif

            execution_plan[i].wave_range_sn = w_range_sn[i,*]   ; Wavelength range for S/N calc
            execution_plan[i].threshold_ston_bin = threshold_ston_bin[i]; Bin S/N limit
            
            execution_plan[i].wave_range_analysis = w_range_analysis[i,*]; Wave range for analysis
            execution_plan[i].threshold_ston_analysis = threshold_ston_analysis[i]; Analysis S/N lim

            ; Set the analysis steps
            execution_plan[i].analysis = $
                    MDAP_SET_EXECUTION_PLAN_ANALYSIS(reform(analysis[i,*]), n_analyses)

            execution_plan[i].analysis_par = analysis_par[i]            ; Analysis parameters

            execution_plan[i].analysis_prior = strcompress(analysis_prior[i],/remove_all)   ; Prior

            ; Check that the prior makes sense
            if strlen(execution_plan[i].analysis_prior) ne 0 then begin
                if file_test(execution_plan[i].analysis_prior) eq 0 then begin
                    indx = fix(execution_plan[i].analysis_prior)

                    ; This is a little contrived because you can't catch
                    ; type conversion errors in IDL, and the returned
                    ; value when fix() fails is 0.  This statement
                    ; checks that the conversion of the fix() value back
                    ; to a string is the same as the input.

                    if MDAP_STC(indx) ne execution_plan[i].analysis_prior then begin
                        message, 'Cannot interpret analysis prior: ' + $
                                 execution_plan[i].analysis_prior
                    endif

;                   if indx lt 0 or indx ge i then begin
                    if indx lt 0 || indx ge i then begin
                        message, 'Plan ' + MDAP_STC(i,/integer) + ' cannot use prior from a plan ' $
                                 + MDAP_STC(indx, /integer) + ' either because there is no such ' $
                                 + 'plan or the ordering is wrong!'
                    endif

                    ; Set the prior to the output file of the designated ExecutionPlan
                    execution_plan[i].analysis_prior = execution_plan[indx].ofile
                endif
            endif

            indx = where(execution_plan[i].analysis eq 1, count)
;           if indx[0] eq -1 then begin; No analysis to perform
            if count eq 0 then begin; No analysis to perform
                execution_plan[i].tpl_lib = -1          ; No template library needed
                execution_plan[i].ems_par = -1          ; No emission-line parameters needed
                execution_plan[i].abs_par = -1          ; No absorption-line parameters needed
            endif else begin
                execution_plan[i].tpl_lib = tpl_lib_analysis[i] ; Template library to use
                execution_plan[i].ems_par = ems_par_analysis[i] ; Emission-line parameters to use
                execution_plan[i].abs_par = abs_par_analysis[i] ; Absorption-line parameters to use
            endelse

            execution_plan[i].overwrite = overwrite_flag[i]     ; Overwrite existing files?

            ; Set the output file name
            execution_plan[i].ofile = MDAP_OUTPUT_FILE_NAME(file_root, execution_plan[i], i)

            ; Perform some (additional) checks of the plan parameters
            execution_plan_i = execution_plan[i]
            MDAP_CHECK_EXECUTION_PLAN, n_tpl, n_ems, n_abs, execution_plan_i
            execution_plan[i] = execution_plan_i
        endfor
END

