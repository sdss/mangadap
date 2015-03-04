;+
; NAME:
;       MDAP_DRP_SNR_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;
;   If the block is to be performed:
;       - Selects good spectra
;       - Determines the median S/N per pixel over the designated
;         wavelength range
;       - Writes the S/N data to the output file
;   Else:
;       - Reads the S/N data to the output file
;
; CALLING SEQUENCE:
;       MDAP_DRP_SNR_BLOCK, manga_dap_version, execution_plan, perform_block, header, wave, flux, $
;                           ivar, mask, gflag, signal, noise, nolog=nolog, $
;                           log_file_unit=log_file_unit, quiet=quiet
;
; INPUTS:
;       manga_dap_version MaNGADAPVersion
;               Structure used to keep track of various
;               version-controlled procedures in the DAP.
;      
;       execution_plan ExecutionPlan
;               Structure providing plan and parameters needed for the
;               analyses.
;
;       perform_block RequiredAnalysisBlock
;               Structure defining which analysis blocks are to be
;               performed.
;
;       header strarr[]
;               Fits header
;            
;       wave dblarr[T]
;               Wavelength of each spectral channel T.
;       
;       flux dblarr[N][T]
;               Galaxy spectra as produced by MDAP_READ_DRP_FITS.
;
;       ivar dblarr[N][T]
;               Inverse variance of the flux
;
;       mask dblarr[N][T]
;               Bad pixel mask.
;
;
; OPTIONAL INPUTS:
;       log_file_unit LUN
;               File unit pointing to the log file
;
; OPTIONAL KEYWORDS:
;       /nolog
;               Suppress output to the log file.
;
;       /quiet
;               Suppress output to the screen.
;
; OUTPUT:
;       gflag intarr[N]
;               Flag (0=false; 1=true) that the spectrum is 'good' as defined by
;               MDAP_SELECT_GOOD_SPECTRA.  Spectra that are NOT good are ignored.
;
;       signal dblarr[N]
;               Mean galaxy signal per angstrom
;
;       noise dblarr[N]
;               Mean galaxy error per angstrom
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
;       MDAP_CALCULATE_SN
;       MDAP_SELECT_GOOD_SPECTRA
;       MDAP_SELECT_WAVE
;       SXADDPAR
;       MDAP_WRITE_OUTPUT
;       MDAP_DEFINE_OUTPUT
;       MDAP_READ_OUTPUT
;
; REVISION HISTORY:
;       01 Feb 2015: Pulled from manga_dap.pro by K. Westfall (KBW)
;       09 Feb 2015: (KBW) Include fraction of good pixels and min(flux)
;                          == max(flux) from MDAP_SELECT_GOOD_SPECTRA in
;                          output file.  Not returned to main.  Theshold
;                          fraction of good pixels required for a "good"
;                          spectrum !! HARD-WIRED !! to 0.8.
;-
;------------------------------------------------------------------------------

PRO MDAP_DRP_SNR_BLOCK, $
        manga_dap_version, execution_plan, perform_block, header, wave, flux, ivar, mask, $
        gflag, signal, noise, nolog=nolog, log_file_unit=log_file_unit, quiet=quiet

        if perform_block.ston eq 1 then begin       ; Perform the S/N calculation if necessary

            ; Get the version
            MDAP_CALCULATE_SN, version=manga_dap_version.calculate_sn

            ; Determine the 'good' spectra based on a set of criteria
            ; defined by this procedure; gflag is returned for later use
            ; (in MDAP_CALCULATE_SN, MDAP_VELOCITY_REGISTER, and
            ; MDAP_SPATIAL_BINNING)
            MDAP_SELECT_GOOD_SPECTRA, flux, ivar, mask, gflag, gindx, good_fraction_threshold=0.8, $
                                      fraction_good=fraction_good, min_eq_max=min_eq_max, $
                                      quiet=quiet

            ; Select the pixels to use in the S/N calculation
            MDAP_SELECT_WAVE, wave, execution_plan.wave_range_sn, lam_sn

            ; Calculate the S/N per pixel over some wavelength range
            MDAP_CALCULATE_SN, flux, ivar, mask, wave, lam_sn, signal, noise, gflag=gflag

            ; Write the version of the S/N calculate ot the header
            SXADDPAR, header, 'VDAPSTON', manga_dap_version.calculate_sn, $
                      'mdap_calculate_sn version'

            ; Add some information to the log
            if ~keyword_set(nolog) then begin
                printf, log_file_unit, '[INFO] S/N calculation over range: ', $
                        execution_plan.wave_range_sn
                printf, log_file_unit, '[INFO] S/N calculation done using version: ' + $
                        manga_dap_version.calculate_sn
            endif

            ; Add the signal and noise to the DRPS extension
            MDAP_WRITE_OUTPUT, execution_plan.ofile, header=header, $
                               w_range_sn=execution_plan.wave_range_sn, signal=signal, $
                               noise=noise, fraction_good=fraction_good, min_eq_max=min_eq_max, $
                               quiet=quiet

        endif else begin                    ; Otherwise read the data from the existing DRP file

            ; Always need the S/N data
            print, 'READING EXISTING S/N DATA'

            MDAP_DEFINE_OUTPUT, header=header, signal=signal, noise=noise
            MDAP_READ_OUTPUT, execution_plan.ofile, header=header, signal=signal, noise=noise

            ; Set gflag based on the values of the signal
            ; TODO: Could now set this using fraction_good and min_eq_max flag instead
            ndrp = n_elements(signal)
            gflag = intarr(ndrp)
            indx = where(signal gt 0, count)
            if count ne 0 then $
                gflag[indx] = 1

        endelse
END

