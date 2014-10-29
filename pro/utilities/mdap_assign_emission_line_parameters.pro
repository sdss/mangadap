;+
; NAME:
;	MDAP_ASSIGN_EMISSION_LINE_PARAMETERS
;
; PURPOSE:
;	Assign arrays previously read in to an array of EmissionLine structures.
;
; CALLING SEQUENCE:
;	eml_par = MDAP_ASSIGN_EMISSION_LINE_PARAMETERS(id, name, lambda, action, kind, a, v, s, fit)
;
; INPUTS:
;	id intarr[E]
;		Unique integer identifyer of the emission line.
;			
;	name strarr[E]
;		Name of the element, read as a string.  Special characters, such
;		as '[', '.' are not permitted.
;			
;	lambda dblarr[E]
;		Rest frame wavelength of the emission line to fit.  WARNING: the
;		name of the emission line in the DAP output is:
;		NAME+'_'+MDAP_STC(ROUND(LAMBDA), /integer)  
;
;	action strarr[E]
;		Describes how the line should be treated.  Possible values are:
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
;	kind strarr[E]
;		Type of line, which can be either 'l' for a line or 'dN' for a
;		doublet.  For doublets, the 'N' indicates the line ID of the
;		other line in the doublet.  The line to which the doublet is
;		tied should have the LINE='l'; for example, if emission line
;		with ID=4 has line='d3', then the emission line with ID=3 must
;		have LINE='l'.
;			
;	A dblarr[E]
;		Relative intensity of the gas emission (positive) or absorption
;		(negative) lines with respect to the doublet.  Therfore, this
;		should most often be unity if KIND='l' and indicate the ratio of
;		line A if KIND='dN'.
;			
;	V dblarr[E]
;		Guess for the velocity offset with respect the galaxy systemic
;		velocity.  TODO: This value is currently ignored by the DAP!
;			
;	S dblarr[E]
;		Guess for the velocity dispersion.  TODO: This value is
;		currently ignored by the DAP!
;			
;	MODE strarr[E]
;		Fitting mode for the line, which can be either 'f' to fit the
;		line independently or 'tN' to set both the velocity and
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
;	30 Sep 2014: (KBW) Copied from MDAP_READ_EMISSION_LINE_PARAMETERS
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_ASSIGN_EMISSION_LINE_PARAMETERS, $
		id, name, lambda, action, kind, a, v, s, fit

	; Declare the array of structures, defining the structure in the process
	eml_par = replicate( { EmissionLine, i:0L, name:'', lambda:0.0d, action:'', kind:'', $
			       a:0.0d, v:0.0d, s:0.0d, fit:'' }, n_elements(id) )

	eml_par.i = id			; Set the IDs
	eml_par.name = name		; Set the Names
	eml_par.lambda = lambda		; Set the Wavelenths
	eml_par.action = action		; Set the Action
	eml_par.kind = kind		; Set the Kind
	eml_par.a = a			; Set the Normalization
	eml_par.v = v			; Set the Velocity
	eml_par.s = s			; Set the Sigma
	eml_par.fit = fit		; Set the Fit mode
	
	return, eml_par
END

