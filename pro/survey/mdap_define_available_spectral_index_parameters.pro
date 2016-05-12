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

; dapsrc is an optional input to define the DAP source path instead of
; using environmental varaibles.

PRO MDAP_DEFINE_AVAILABLE_SPECTRAL_INDEX_PARAMETERS, $
        abs_line_keys, absorption_line_parameters, abs_vacuum_wave, dapsrc=dapsrc

        ; Define the DAP source path
        if n_elements(dapsrc) eq 0 then $
            dapsrc = getenv('MANGADAP_DIR')

        ;-----------------------------------------------------------------------
        ; Define the set of absorption-line parameter files.  The format expected
        ; for these files is described above.
        nabs_files = 1
        abs_line_keys = strarr(nabs_files)
        absorption_line_parameters = strarr(nabs_files)
        abs_vacuum_wave = intarr(nabs_files)

        abs_line_keys[0] = 'LICK'
        absorption_line_parameters[0] = dapsrc+'/external/legacy/absorption_line_indices_definition.dat'
        abs_vacuum_wave[0] = 0

END

