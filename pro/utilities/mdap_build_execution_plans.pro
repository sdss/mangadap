;+
; NAME:
;       MDAP_BUILD_EXECUTION_PLANS
;
; PURPOSE:
;       Populate a set of ExecutionPlan structures defining how the DAP is
;       executed for a given DRP-produced fits file.
;
; CALLING SEQUENCE:
;       MDAP_BUILD_EXECUTION_PLANS, n_tpl, n_ems, n_abs, bin_type, w_range_sn, bin_par, $
;                                   threshold_ston_bin, bin_weight_by_sn2, w_range_analysis, $
;                                   threshold_ston_analysis, analysis, tpl_lib_analysis, $
;                                   ems_par_analysis, abs_par_analysis, overwrite_flag, $
;                                   execution_plan
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
;       bin_type strarr[P]
;               Type of binning to apply for each of the P plans to create.
;               Valid values are:
;                       'NONE' - No binning is performed, all spectra are
;                                analyzed.
;                       'ALL' - All spectra are combined into a single bin for
;                               analysis.
;                       'STON' - Spectra are binned, using the Voronoi binning
;                                scheme, to a minimum S/N.  This type require a
;                                single parameter (provided via bin_par), which
;                                is the minimum S/N level.
;                       'RAD' - Spectra are binned in radius.  This type
;                               requires five parameters, the x and y center,
;                               ellipticity, position angle, and width of the
;                               radial bins.  
;                       
;                               TODO: When binning, the spectra are deredshifted
;                               before creating the bin, meaning that a set of
;                               velocities for each spectrum must be available!
;                       
;                               TODO: Currently only the bin width is input!!
;                       
;       w_range_sn dblarr[P][2]
;               Wavelength range used to calculate the S/N of each spectrum,
;               performed using each individual spectrum.  This is always
;               calculated, regardless of the binning type.
;
;       bin_par dblarr[P]
;               Parameters required for the binning.
;               TODO: Can I make this a structure?
;
;       threshold_ston_bin dblarr[P]
;               The S/N threshold for the inclusion of any DRP spectrum in the
;               binning process.
;
;       bin_weight_by_sn2 intarr[P]
;               Flag to weight each spectrum by S/N^2 when producing the binned
;               spectra.  If true (eq 1), weights are S/N^2; if false (eq 0),
;               the spectra are added with uniform weights.
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
;                       'stellar-cont' - Fit the stellar continuum.  Requires
;                                        the template library.  Will mask
;                                        emission lines, if a set of
;                                        emission-line parameters are provided.
;                                        Determines the optimal template mix and
;                                        the stellar kinematics.
;
;                       'emission-line' - Fits the emission lines.  Requires a
;                                         template library and a set of
;                                         emission-line parameters and a
;                                         previous fit to the stellar
;                                         kinematics.  Will re-optimize the
;                                         template mix, but does not adjust the
;                                         stellar kinematics.  Determines the
;                                         emission-line kinematics, intensities,
;                                         fluxes, and equivalent widths.
;
;                       'abs-indices' - Calculates the absorption-line indices.
;                                       Requires an absorption-line parameter
;                                       set and a previous fit to the stellar
;                                       continuum.  TODO: Also requires fit to
;                                       emission-line parameters?
;
;       analysis_extra strarr[P][]
;               List of extra parameters to set for each of the P analyses to
;               perform.  Zero length strings are ignored.
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
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       execution_plan ExecutionPlan[P]
;               Array of P structures containing the execution plans.
;
; OPTIONAL OUTPUT:
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
;       10 Oct 2014: (KBW) Original implementation
;       13 Oct 2014: (KBW) Include bin_weight_by_sn2
;       15 Oct 2014: (KBW) Added analysis_extra
;       07 Nov 2014: (KBW) Checks that analysis steps are performed in
;                          sequential order is now done using
;                          MDAP_ANALYSIS_BLOCKS_TO_PERFORM(), which takes into
;                          account the blocks that have already been completed.
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_SET_EXECUTION_PLAN_ANALYSIS, $
                analysis

        n_steps = n_elements(analysis)          ; Number of analysis steps

        analysis_ = intarr(3)                   ; Currently only 3 possible analysis steps
        if n_steps eq 0 then $
            return, analysis_                   ; All analysis steps set to zero

        ; Set the steps
        for i=0,n_steps-1 do begin
            if analysis[i] eq 'stellar-cont' then begin
                analysis_[0] = 1                ; Perform the stellar continuum fitting
                continue
            endif
            if analysis[i] eq 'emission-line' then begin
                analysis_[1] = 1                ; Perform the emission line fitting
                continue
            endif
            if analysis[i] eq 'abs-indices' then begin
                analysis_[2] = 1                ; Perform the absorption index measurements
                continue
            endif
            
            if strlen(analysis[i]) ne 0 then $
                message, 'Unknown analysis step: ', analysis[i]
        endfor

        ; TODO: So far all steps must occur serially.

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
        if execution_plan.bin_type ne 'NONE' and execution_plan.bin_type ne 'ALL' and $
           execution_plan.bin_type ne 'STON' and execution_plan.bin_type ne 'RAD' then begin 
            message, 'Unknown binning type:' + execution_plan.bin_type
        endif

        ; Check bin parameters
        if (execution_plan.bin_type eq 'STON' and n_elements(execution_plan.bin_par) ne 1) or $
           (execution_plan.bin_type eq 'RAD' and n_elements(execution_plan.bin_par) ne 1) then begin
            message, 'STON and RAD bin types must have one binning parameter (S/N and Bin width,'+ $
                   ' respectively)'
            
        endif

        ; Analysis steps are check using MDAP_SET_EXECUTION_PLAN_ANALYSIS

        ; Check that the selected template library exists
        if execution_plan.tpl_lib ne -1 and execution_plan.tpl_lib ge n_tpl then $
            message, 'No template library index: ', execution_plan.tpl_lib
            
        ; Check that the selected emission-line parameter set exists
        if execution_plan.ems_par ne -1 and execution_plan.ems_par ge n_ems then $
            message, 'No emission-line parameter set: ', execution_plan.ems_par
            
        ; Check that the selected absorption-line parameter set exists
        if execution_plan.abs_par ne -1 and execution_plan.abs_par ge n_abs then $
            message, 'No absorption-line parameter set: ', execution_plan.abs_par


        ; If analysis parameters are not needed, then turn them off
        ;
        ; Absorption-line parameters only needed for that analysis
        if execution_plan.analysis[2] eq 0 then $
            execution_plan.abs_par = -1
        ; Emission-line parameters used for both stellar-continuum and emission-line analysis
        if execution_plan.analysis[1] eq 0 and execution_plan.analysis[0] eq 0 then $
            execution_plan.ems_par = -1
        ; Template spectra used for absorption-line, emission-line, and stellar continuum analysis
        if execution_plan.analysis[2] eq 0 and execution_plan.analysis[1] eq 0 and $
            execution_plan.analysis[0] eq 0 then begin
            execution_plan.tpl_lib = -1
        endif
END

PRO MDAP_BUILD_EXECUTION_PLANS, $
                n_tpl, n_ems, n_abs, bin_type, w_range_sn, bin_par, threshold_ston_bin, $
                bin_weight_by_sn2, w_range_analysis, threshold_ston_analysis, analysis, $
                analysis_par, tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                overwrite_flag, execution_plan

        n_plans = n_elements(bin_type)                  ; Number of plans to execute
        n_extra = (size(analysis_extra))[2]             ; Maximum number of analysis extras

        analysis_par_def = MDAP_DEFINE_ANALYSIS_PAR()

        ; Instantiate the ExecutionPlan structure
        ; TODO: Make analysis_extra a pointer?
        execution_plan = replicate( { ExecutionPlan, bin_type:'', wave_range_sn:dblarr(2), $
                                                     bin_par:0.0d, threshold_ston_bin:0.0d, $
                                                     bin_weight_by_sn2:0, $
                                                     wave_range_analysis:dblarr(2), $
                                                     threshold_ston_analysis:0.0d, $
                                                     analysis:intarr(3), $
                                                     analysis_par:analysis_par_def, tpl_lib:0, $
                                                     ems_par:0, abs_par:0, overwrite:0, $
                                                     ofile:''}, n_plans)

        for i=0,n_plans-1 do begin

            execution_plan[i].bin_type = bin_type[i]            ; Bin type
            execution_plan[i].wave_range_sn = w_range_sn[i,*]   ; Wavelength range for S/N calc

            execution_plan[i].bin_par = bin_par[i]                      ; Binning parameters
            execution_plan[i].threshold_ston_bin = threshold_ston_bin[i]; Bin S/N limit
            execution_plan[i].bin_weight_by_sn2 = bin_weight_by_sn2[i]  ; Flag to use S/N^2 weights
            
            execution_plan[i].wave_range_analysis = w_range_analysis[i,*]; Wave range for analysis
            execution_plan[i].threshold_ston_analysis = threshold_ston_analysis[i]; Analysis S/N lim

            ; Set the analysis steps
            execution_plan[i].analysis = MDAP_SET_EXECUTION_PLAN_ANALYSIS(reform(analysis[i,*]))

            execution_plan[i].analysis_par = analysis_par[i]    ; Analysis parameters

            indx = where(execution_plan[i].analysis eq 1)
            if indx[0] eq -1 then begin; No analysis to perform
                print, 'bad index'
                execution_plan[i].tpl_lib = -1          ; No template library needed
                execution_plan[i].ems_par = -1          ; No emission-line parameters needed
                execution_plan[i].abs_par = -1          ; No absorption-line parameters needed
            endif else begin
                execution_plan[i].tpl_lib = tpl_lib_analysis[i] ; Template library to use
                execution_plan[i].ems_par = ems_par_analysis[i] ; Emission-line parameters to use
                execution_plan[i].abs_par = abs_par_analysis[i] ; Absorption-line parameters to use
            endelse

            execution_plan[i].overwrite = overwrite_flag[i]     ; Overwrite existing files?

            execution_plan_i = execution_plan[i]
            MDAP_CHECK_EXECUTION_PLAN, n_tpl, n_ems, n_abs, execution_plan_i
            execution_plan[i] = execution_plan_i
        endfor
END

