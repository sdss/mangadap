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

; dapsrc is an optional input to define the DAP source path instead of
; using environmental varaibles.

PRO MDAP_DEFINE_AVAILABLE_EMISSION_LINE_PARAMETERS, $
        ems_line_keys, emission_line_parameters, ems_vacuum_wave, dapsrc=dapsrc

        ; Define the DAP source path
        if n_elements(dapsrc) eq 0 then $
            dapsrc = getenv('MANGADAP_DIR')

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

END

