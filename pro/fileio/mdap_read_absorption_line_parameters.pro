;+
; NAME:
;	MDAP_READ_ABSORPTION_LINE_PARAMETERS
;
; PURPOSE:
;
; CALLING SEQUENCE:
;	indices = MDAP_READ_ABSORPTION_LINE_PARAMETERS(file)
;
; INPUTS:
;	file string
;		Name of the file containing the indices to measure.  File has 9
;		columns:
;
;	# ID Name        pass_0    pass_1   blcnt_0   blcnt_1   rdcnt_0   rdcnt_1 units
;	#                   ang       ang       ang       ang       ang       ang
;  	   0 D4000      0000.00   0000.00   3750.00   3950.00   4050.00   4250.00   ang 
;	   1 CaII0p39   3899.50   4003.50   3806.50   3833.80   4020.70   4052.40   ang
;	   2 HDeltaA    4083.50   4122.25   4041.60   4079.75   4128.50   4161.00   ang
;	   3 HDeltaF    4091.00   4112.25   4057.25   4088.50   4114.75   4137.25   ang
;
;		1 - ID: ID number for the index (integer)
;
;		2 - Name: Name for the index (string)
;
;		3,4 - pass_0,pass_1: Starting and ending wavelength over which
;		to get the equivalent width (TODO: CHECK). (float)
;		
;		5,6 - blcnt_0,blcnt_1: Starting and ending wavelength over which
;		to measure the continuum level blueward of the line index. (float)
;		
;		7,8 - rdcnt_0,rdcnt_1: Starting and ending wavelength over which
;		to measure the continuum level redward of the line index. (float)
;		
;		9 - units: Units to return for the line index
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	indices AbsorptionIndex[]
;		An array of structures that contain the defining parameters of
;		each absorption-line index.
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
;	MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	03 Nov 2014: (KBW) Copied from L. Coccato (v0_8)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_READ_ABSORPTION_LINE_PARAMETERS, $
		file

	; Read the data
	READCOL, file, id, name, pass0, pass1, blcnt0, blcnt1, rdcnt0, rdcnt1, units, $
		format='I,A,D,D,D,D,D,D,A', comment='#', /silent

	; Assign the parameters
	return, MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS(id, name, pass0, pass1, blcnt0, blcnt1, $
						       rdcnt0, rdcnt1, units)
END


