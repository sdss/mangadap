;===============================================================================
; STELLAR TEMPLATE LIBRARY FILES
;===============================================================================
;	The stellar template library files should be a list of 1D fits files,
;	and be associated with one of the following library keys (TODO: More to
;	be added):
;
;	                         Spectral
;	          KEY    resolution (ang)	
;	-------------    ----------------
;	    M11-MARCS                2.73
;
;	The fits files must have CRVAL1, CRPIX1, and CDELT1 keywords used to
;	define the wavelength coordinates of each pixel:
;
;		pix = (findgen(npixels+1))[1:*]
;		wave = (pix - CRPIX1) * CDELT1 + CRVAL1
;
;	The reference frame of the template wavelengths must also be defined as
;	either vacuum or air via the tpl_vacuum_wave vector.  Set the value to
;	be 1(true) if the wavelengths are in vacuum, 0(false) otherwise.  It is
;	expected that the DRP spectra are in vacuum wavelengths.  The DAP will
;	therefore use the IDL routine AIRTOVAC to convert the template
;	wavelengths to vacuum if tpl_vacuum_wave = 0.
;
;===============================================================================

;===============================================================================
; EMISSION-LINE PARAMETER FILES
;===============================================================================
;		File name with information regarding the emission lines to fit.
;		The format must be compatible with the module that performs the
;		fit of the emission lines.  If Gandalf or S-Gandalf are used,
;		the input file must be an ascii file with 9 columns as in the
;		following example (comments starts with ''\#''):
;
;	#  ID     CODE    WAVELENGTH   ACTION    LINE/  INTENSITY   V_g/i   sig_g/i   MODE
;	#                 [angstrom]  i/f/m/s  doublet                                f/tN
;	#   0     HeII       3203.15        m        l      1.000       0       10     t25
;	#   1     NeV        3345.81        m        l      1.000       0       10     t25
;	    2     NeV        3425.81        m        l      1.000       0       10     t25
;	    3     OII        3726.03        m        l      1.000       0       10     t25
;              
;		The columns are:
;
;		1. ID: Unique integer identifyer of the emission line.
;			
;		2. CODE: Name of the element, read as a string.  Special
;		characters, such as '[', '.' are not permitted.
;			
;		3. WAVELENGTH: Rest frame wavelength of the emission line to
;		fit.  WARNING: the name of the emission line in the DAP output
;		is: CODE+'_'+MDAP_STC(ROUND(wav), /integer)  
;			
;		4. ACTION: Describes how the line should be treated.  Possible
;		values are:
;
;			'i': ignore the line, as if the line were commented out.
;
;			'f': fit the line and mask the line when fitting the
;			stellar continuum.
;
;			'm': mask the line when fitting the stellar continuum
;			but do NOT fit the line itself
;
;			's': defines a sky line that should be masked.  When
;			masked, the wavelength of the line is NOT adjusted for
;			the redshift of the object spectrum.
;			
;		5. LINE:  Type of line, which can be either 'l' for a line or
;		'dN' for a doublet.  For doublets, the 'N' indicates the line ID
;		of the other line in the doublet.  The line to which the doublet
;		is tied should have the LINE='l'; for example, if emission line
;		with ID=4 has line='d3', then the emission line with ID=3 must
;		have LINE='l'.
;			
;		6. INTENSITY:  Relative intensity of the gas emission (positive)
;		or absorption (negative) lines with respect to the doublet.
;		Therfore, this should most often be unity if LINE='l' and
;		indicate the ratio of line INTENSITY if LINE='dN'.
;			
;		7. V_g/i: Guess for the velocity offset with respect the galaxy
;		systemic velocity.  TODO: This value is currently ignored by the
;		DAP!
;			
;		8. sig_g/i. Guess for the velocity dispersion.  TODO: This value
;		is currently ignored by the DAP!
;			
;		9. MODE.  Fitting mode for the line, which can be either 'f' to
;		fit the line independently or 'tN' to set both the velocity and
;		dispersion to be tied to a fitted line (MODE='f') with ID=N.
;		One can also force only the velocities to be tied using 'vN' or
;		only the velocity dispersions to be tied using 'sN'.
;
;	The reference frame of the emission-line wavelengths must also be
;	defined as either vacuum or air via the ems_vacuum_wave vector.  Set the
;	value to be 1(true) if the wavelengths are in vacuum, 0(false)
;	otherwise.  It is expected that the DRP spectra are in vacuum
;	wavelengths.  The DAP will therefore use the IDL routine AIRTOVAC to
;	convert the emission-line wavelengths to vacuum if ems_vacuum_wave = 0.
;
;===============================================================================

;===============================================================================
;===============================================================================

;===============================================================================
; SPECTRAL INDEX PARAMETER FILES
;===============================================================================
;	Files must be ascii text with the following 9 columns, e.g.:
;
;	# ID Name        pass_0    pass_1   blcnt_0   blcnt_1   rdcnt_0   rdcnt_1 units
;	#                   ang       ang       ang       ang       ang       ang
;  	   0 D4000      0000.00   0000.00   3750.00   3950.00   4050.00   4250.00   ang 
;	   1 CaII0p39   3899.50   4003.50   3806.50   3833.80   4020.70   4052.40   ang
;	   2 HDeltaA    4083.50   4122.25   4041.60   4079.75   4128.50   4161.00   ang
;	   3 HDeltaF    4091.00   4112.25   4057.25   4088.50   4114.75   4137.25   ang
;
;		1. Integer. Unique ID number of the absorption line feature.
;
;		2. String. Unique name of the absorption line feature. This will
;		define the name of the field in sctructure of the DAP results
;		(i.e. the name must begin with a letter, special characters like
;		comas or dots are not allowed).
;
;		3-4. Float (units: ang) Lower and upper value of the index
;		passband.  If these two number are both less than one (see
;		MDAP_INDEX_IS_BANDHEAD), the index is treated as a bandhead or
;		spectral break.
;
;		5-6. Float (units: ang) Lower and upper value of the index blue
;		pseudo-continuum.
;
;		7-8. Float (units: ang) Lower and upper value of the index red
;		pseudo-continuum.
;
;		9. String (accepted values are: ang or mag). Specifies the units
;		(ang or magnitudes) of the output.
;
;	The reference frame of the absorption-line wavelength parameters must
;	also be defined as either vacuum or air via the abs_vacuum_wave vector.
;	Set the value to be 1(true) if the wavelengths are in vacuum, 0(false)
;	otherwise.  It is expected that the DRP spectra are in vacuum
;	wavelengths.  The DAP will therefore use the IDL routine AIRTOVAC to
;	convert the absorption-line wavelength parameters to vacuum if
;	ems_vacuum_wave = 0.
;
;	NOTE: Indices will be measured only if their blue and red
;	pesudo-continua bandpasses are included in the considered wavelength
;	range. If not, their values are set to NaN, and their errors to 99 in
;	the final output file.
;
;===============================================================================
;===============================================================================

;===============================================================================
; SPATIAL BINNING TYPES
;===============================================================================
;	Type of binning to apply for each of the P plans to create.
;	Valid values are:
;		'NONE' - No binning is performed, all spectra are
;		         analyzed.
;		'ALL' - All spectra are combined into a single bin for
;			analysis.
;		'STON' - Spectra are binned, using the Voronoi binning
;			 scheme, to a minimum S/N.  This type require a
;			 single parameter (provided via bin_par), which
;			 is the minimum S/N level.
;
;	TODO: Need to add 'RAD' option for radial binning
;===============================================================================
;===============================================================================

;===============================================================================
; ANALYSIS STEPS
;===============================================================================
;	The currently available analysis steps are:
;
;		'stellar-cont' - Fit the stellar continuum.  Requires
;				 the template library.  Will mask
;				 emission lines, if a set of
;				 emission-line parameters are provided.
;				 Determines the optimal template mix and
;				 the stellar kinematics.
;
;		'emission-line' - Fits the emission lines.  Requires a
;				  template library and a set of
;				  emission-line parameters and a
;				  previous fit to the stellar
;				  kinematics.  Will re-optimize the
;				  template mix, but does not adjust the
;				  stellar kinematics.  Determines the
;				  emission-line kinematics, intensities,
;				  fluxes, and equivalent widths.
;
;		'abs-indices' - Calculates the absorption-line indices.
;				Requires an absorption-line parameter
;				set and a previous fit to the stellar
;				continuum.  TODO: Also requires fit to
;				emission-line parameters?
;
;				TODO: Absorption-line analysis currently not
;				performed!
;
; NOTE: The current list of fits must occur sequentially.  That is, you cannot
; pick only 'emission-line' analysis.  If you do, the program will warn you and
; then also perform the 'stellar-cont' analysis.
;===============================================================================

; Setup some necessary execution variables for the MANGA DAP

; Signifier is what will be reported to the log file as the configuration file
; for a run of MaNGA_DAP

PRO MDAP_EXECUTION_SETUP, $
	total_filelist, output_root_dir, tpl_library_keys, template_libraries, tpl_vacuum_wave, $
	ems_line_keys, emission_line_parameters, ems_vacuum_wave, abs_line_keys, $
	absorption_line_parameters, abs_vacuum_wave, signifier, bin_type, bin_par, w_range_sn, $
	threshold_ston_bin, bin_weight_by_sn2, w_range_analysis, threshold_ston_analysis, $
	analysis, tpl_lib_analysis, ems_par_analysis, abs_par_analysis, overwrite_flag, $
	analysis_par, save_intermediate_steps=save_intermediate_steps, $
	remove_null_templates=remove_null_templates, external_library=external_library

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
	
	total_filelist = './mdap_all_inp_nsa_v1b_0_0_v2.par' 

	;-----------------------------------------------------------------------
	; Location for the output files
	;
	;	For a file name of 
	;		manga-7443-12701-LOGCUBE.fits
	;
	;	Output files will be placed in
	;		output_root_dir+'/'+manga-7443-12701-LOGCUBE

	output_root_dir=getenv('MANGA_SPECTRO_ANALYSIS')

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

	; external_library=getenv('MANGA_DAP_DIR')+'/external/F90_32/'
	external_library=getenv('MANGA_DAP_DIR')+'/external/F90_64/'

	;-----------------------------------------------------------------------
	; Define the set of template libraries.  The format expected for these
	; files is described above.
	ntpl_libraries = 1
	tpl_library_keys = strarr(ntpl_libraries)
	template_libraries = strarr(ntpl_libraries)
	tpl_vacuum_wave = intarr(ntpl_libraries)

	tpl_library_keys[0] = 'M11-MARCS'
	template_libraries[0] = getenv('MANGA_DAP_DIR')+'/external/templates_marcs/*_s.fits'
	tpl_vacuum_wave[0] = 0

	;-----------------------------------------------------------------------
	; Define the set of emission-line parameter files.  The format expected
	; for these files is described above.
	neml_files = 3
	ems_line_keys = strarr(neml_files)
	emission_line_parameters = strarr(neml_files)
	ems_vacuum_wave = intarr(neml_files)

	ems_line_keys[0] = 'STANDARD'
	emission_line_parameters[0] = getenv('MANGA_DAP_DIR') + $
			'/external/emission_lines_setup_with_Balmer_decrement'
	ems_vacuum_wave[0] = 0

	ems_line_keys[1] = 'NODOUBLETS'
	emission_line_parameters[1] = getenv('MANGA_DAP_DIR') + $
			'/external/emission_lines_setup_with_Balmer_decrement_no_doublets'
	ems_vacuum_wave[1] = 0

	ems_line_keys[2] = 'RESIDUAL'
	emission_line_parameters[2] = getenv('MANGA_DAP_DIR') + $
			'/external/emission_lines_setup_with_Balmer_decrement_residuals'
	ems_vacuum_wave[2] = 0

	;-----------------------------------------------------------------------
	; Define the set of absorption-line parameter files.  The format expected
	; for these files is described above.
	nabs_files = 1
	abs_line_keys = strarr(nabs_files)
	absorption_line_parameters = strarr(nabs_files)
	abs_vacuum_wave = intarr(nabs_files)

	abs_line_keys[0] = 'LICK'
	absorption_line_parameters[0] = getenv('MANGA_DAP_DIR') + $
			'/external/absorption_line_indices_definition.dat'
	abs_vacuum_wave[0] = 0

	;=======================================================================
	; DEFINITION OF EXECUTION PROCEDURES

	; Define a string used to signify this file in the header of the DAP output file(s)
	cd, current=directory
	signifier = directory+'/mdap_setup.pro'

	;-----------------------------------------------------------------------
	; Define the number of execution iterations and setup the needed vectors
	; and allocate the necessary arrays.

	niter = 1					; Number of ExecutionPlans to produce

	bin_type = strarr(niter)			; Binning type
	bin_par = dblarr(niter)				; Binning parameter

	w_range_sn = dblarr(niter, 2)			; Wavelength range for S/N calculation
	threshold_ston_bin = dblarr(niter)		; Threshold S/N to include spectrum in bin
	bin_weight_by_sn2 = intarr(niter)		; Weight by S/N^2 when combining spectra

	w_range_analysis = dblarr(niter, 2)		; Wavelength range for the analysis
	threshold_ston_analysis = dblarr(niter)		; Threshold S/N to analyze spectrum
	max_analysis_blocks = 3				; Maximum number of analysis blocks
	analysis = strarr(niter, max_analysis_blocks)	; Analysis steps to apply

	tpl_lib_analysis = intarr(niter)		; INDEX of template library to use
	ems_par_analysis = intarr(niter)		; INDEX of emission-line parameter file
	abs_par_analysis = intarr(niter)		; INDEX of absorption-line parameter file

	; A structure to hold parameters required for the analysis
	analysis_par_def = MDAP_DEFINE_ANALYSIS_PAR()
	analysis_par = replicate( analysis_par_def, niter)

	overwrite_flag = intarr(niter)			; Flag to overwrite any existing output file

	;-----------------------------------------------------------------------
	; For each iteration:
	
	; Define the binning type.  The available binning types are listed
	; above.
	bin_type[0] = 'STON'

	; Define the bin parameter. See the discussion of the binning types above
	bin_par[0] = 40.0d

	; Define the wavelength range over which to calculate the mean S/N per pixel
	w_range_sn[0,*] = [5560.00, 6942.00]

	; Define the S/N threshold to include spectrum in any bin
	threshold_ston_bin[0] = -300.0d

	; Set the flag for whether or not to bin by S/N^2 (1-yes, 0-no)
	bin_weight_by_sn2[0] = 1

	; Define the wavelength range over which to perform ALL analyses
	w_range_analysis[0,*] = [3650.,10300.] 

	; Define the S/N threshold to perform analysis
	threshold_ston_analysis[0] = 0.0d

	; Set the list of analyses to perform.  The available analysis steps are
	; listed above.

	analysis[0,0] = 'stellar-cont'
	analysis[0,1] = 'emission-line'
	analysis[0,2] = 'abs-indices'

;	analysis[0,0] = 'stellar-cont'
;	analysis[0,1] = 'abs-indices'

;	analysis[0,0] = 'abs-indices'

	; Set the index of the template library to use for the analysis
	; TODO: Change this to use the library key?
	tpl_lib_analysis[0] = 0

	; Set the index of the emission-line parameter set to use
	ems_par_analysis[0] = 0

	; Set the index of the absorption-line parameter set to use
	abs_par_analysis[0] = 0

	; Set additional parameters needed by the analysis modules
	; The reddening order can be 0, 1, or 2
	; TODO: Allow for oversample?
	; IF NOT SET HERE, the default values are:
	;	moments=2, degree=-1, mdegree=-1, reddening_order=0
	analysis_par[0].moments = 4
	analysis_par[0].degree = -1
	analysis_par[0].mdegree = 6
	analysis_par[0].reddening_order = 0
	; analysis_par[0].reddening[*] = [0.01,0.01]

	; Set a flag to overwrite existing output: 1-yes, 0-no
	overwrite_flag[0] = 0

	;=======================================================================
	;=======================================================================

END

