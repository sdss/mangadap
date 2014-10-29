;+
; NAME:
;	MDAP_READ_EMISSION_LINE_PARAMETERS
;
; PURPOSE:
;	Read the necessary emission line parameters from a file and put the
;	result into a structure.
;
; CALLING SEQUENCE:
;	eml_par = MDAP_READ_EMISSION_LINE_PARAMETERS(file)
;
; INPUTS:
;	file string
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
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	eml_par EmissionLine[E]
;		The parameters for each of E emission lines.  The EmissionLine
;		structure is defined as follows:
;
;		{ EmissionLine, i:0L, name:'', lambda:0.0d, action:'', $
;		      kind:'', a:0.0d, v:0.0d, s:0.0d, fit:'' }
;
;		Once created, one selects, for example, the name of line i=2
;		using: eml_par[2].name
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- This is NOT the same as used in GANDALF, but it may actually function
;	  the same.  I prefer to define a structure that defines one line, than
;	  use one structure that defines all lines.  The former is the
;	  object-oriented way of thinking about things.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	25 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_READ_EMISSION_LINE_PARAMETERS, $
		file

	; Check that the file exists
	if file_test(file) eq 0 then $
	    message, file+' does not exist!'

	; Read the data
	READCOL, file, id, name, lambda, action, kind, a, v, s, fit, FORMAT='(I,A,F,A,A,F,F,F,A)', $
		 COMMENT='#', /SILENT

	return, MDAP_ASSIGN_EMISSION_LINE_PARAMETERS(id, name, lambda, action, kind, a, v, s, fit)
END

