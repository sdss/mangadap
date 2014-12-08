;+
; NAME:
;       MDAP_PRINT_EXECUTION_PLAN
;
; PURPOSE:
;       Report an execution plan to the screen.
;
; CALLING SEQUENCE:
;       MDAP_PRINT_EXECUTION_PLAN, template_libraries, emission_line_parameters, 
;                                  absorption_line_parameters, execution_plan
;
; INPUTS:
;       template_libraries strarr[]
;               The list of available template libraries.
;
;       emission_line_parameters strarr[]
;               The list of available emission-line parameter files.
;
;       absorption_line_parameters strarr[]
;               The list of available absorption-line parameter files.
;
;       execution_plan ExecutionPlan[P]
;               An array of P structures containing an execution plan.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
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
;       13 Oct 2014: (KBW) Include flag for weighting by S/N^2
;       04 Dec 2014: (KBW) Edits to print the now structure
;                          execution_plan.bin_par.
;       05 Dec 2014: (KBW) Print the prior
;-
;------------------------------------------------------------------------------

PRO MDAP_PRINT_EXECUTION_PLAN, $
                template_libraries, emission_line_parameters, absorption_line_parameters, $
                execution_plan
        print, '    Binning type: ', execution_plan.bin_par.type
        if execution_plan.bin_par.type ne 'NONE' then begin
            print, '    Velocity register the spectra before binning: ', $
                   execution_plan.bin_par.v_register
            print, '    Use S/(N)^2 weighting when combining spectra: ', $
                   execution_plan.bin_par.optimal_weighting
        endif
        if execution_plan.bin_par.type eq 'STON' then $
            print, '    Minimum S/N per bin: ', execution_plan.bin_par.ston
        if execution_plan.bin_par.type eq 'RADIAL' then begin
            print, '    X center: ', execution_plan.bin_par.cx
            print, '    Y center: ', execution_plan.bin_par.cy
            print, '    Position angle: ', execution_plan.bin_par.pa
            print, '    Ellipticity: ', execution_plan.bin_par.ell
            print, '    Starting radius: ', execution_plan.bin_par.rs
            print, '    Ending radius (-1 for maximum): ', execution_plan.bin_par.re
            print, '    Number of radial bins: ', execution_plan.bin_par.nr
            print, '    Logarithmic bins: ', execution_plan.bin_par.rlog
            print, '    Radius scale parameter: ', execution_plan.bin_par.rscale
        endif

        print, '    Wavelength range for S/N calculation: ', execution_plan.wave_range_sn
        print, '    S/N threshold for inclusion in bin: ', execution_plan.threshold_ston_bin
        print, '    Wavelength range for analysis: ', execution_plan.wave_range_analysis
        print, '    S/N threshold for inclusion in analysis: ', $
               execution_plan.threshold_ston_analysis

        indx=where(execution_plan.analysis eq 1)
        if indx[0] eq -1 then begin
            print, '    Analyses to complete: NONE'
        endif else begin
            print, '    Analyses to complete: '
            if execution_plan.analysis[0] eq 1 then $
                print, '                          stellar-cont'
            if execution_plan.analysis[1] eq 1 then $
                print, '                          emission-line'
            if execution_plan.analysis[2] eq 1 then $
                print, '                          abs-indices'

            print, '    Analysis parameters: '
            print, '              moments: ', execution_plan.analysis_par.moments
            print, '               degree: ', execution_plan.analysis_par.degree
            print, '              mdegree: ', execution_plan.analysis_par.mdegree
            print, '      reddening order: ', execution_plan.analysis_par.reddening_order
            if execution_plan.analysis_par.reddening_order gt 0 then begin $
                str='[' +MDAP_STC(execution_plan.analysis_par.reddening[0])
                if execution_plan.analysis_par.reddening_order eq 2 then $
                    str=str+','+MDAP_STC(execution_plan.analysis_par.reddening[1])
                str=str+']'
                print, '      reddening coeff: ', str
            endif

            if strlen(execution_plan.analysis_prior) ne 0 then $
                print, '    Analysis Prior: '+execution_plan.analysis_prior

;           n_extras = n_elements(execution_plan.analysis_extra)
;           for i=0,n_extras-1 do begin
;               if strlen(execution_plan.analysis_extra[i]) gt 0 then begin
;                   print, 'Analysis extra ' + strcompress(string(i+1),/remove_all) + ': ' + $
;                          execution_plan.analysis_extra[i]
;               endif
;           endfor

            if execution_plan.tpl_lib ne -1 then $
                print, '    Template library: ', template_libraries[execution_plan.tpl_lib]
            if execution_plan.ems_par ne -1 then begin
                print, '    Emission-line parameters: ', $
                        emission_line_parameters[execution_plan.ems_par]
            endif
            if execution_plan.abs_par ne -1 then begin
                print, '    Absorption-line parameters: ', $
                        absorption_line_parameters[execution_plan.abs_par]
            endif
        endelse
        print, '    Overwrite flag: ', execution_plan.overwrite
        print, '    Output file: ', execution_plan.ofile
END


