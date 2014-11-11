;+
; NAME:
;	MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS
;
; PURPOSE:
;	Create and assign the values for a set of spectral index
;	(absorption-line) parameters.
;
; CALLING SEQUENCE:
;	result = MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS(id, name, pass0, pass1, blcnt0, blcnt1, $
;							rdcnt0, rdcnt1, units)
;
; INPUTS:
;	id intarr[I]
;		Identifying number for each of the I spectral indices.
;
;	name strarr[I]
;		String representing a unique name for the I spectral indices.
;
;	pass0 dblarr[I]
;		Starting wavelength of the main passband of the spectral index.
;		Should be 0.0 and equal to the end of the passband for
;		measurements of a spectrum break or bandhead.
;
;	pass1 dblarr[I]
;		Ending wavelength of the main passband of the spectral index.
;		Should be 0.0 and equal to the end of the passband for
;		measurements of a spectrum break or bandhead.
;
;	blcnt0 dblarr[I]
;		Starting wavelength of the spectral region blueward of the index
;		to use for the pseudo continuum measurement.
;
;	blcnt1 dblarr[I]
;		Ending wavelength of the spectral region blueward of the index
;		to use for the pseudo continuum measurement.
;
;	rdcnt0 dblarr[I]
;		Starting wavelength of the spectral region redward of the index
;		to use for the pseudo continuum measurement.
;
;	rdcnt1 dblarr[I]
;		Ending wavelength of the spectral region redward of the index to
;		use for the pseudo continuum measurement.
;
;	units strarr[I]
;		String representation of the units of the index measurement.
;		Should be 'ang' or 'mag'.
;
; OUTPUT:
;	The result is an array of I SpectralIndex structures with the input
;	properties.
;
; REVISION HISTORY:
;	06 Nov 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS, $
		id, name, pass0, pass1, blcnt0, blcnt1, rdcnt0, rdcnt1, units

	nindex = n_elements(id)					; Number of index definitions

	; Define the structure
	indices = replicate( { SpectralIndex, id:0, name:'', units:'', passband:dblarr(2), $
					      blue_cont:dblarr(2), red_cont:dblarr(2) }, nindex)

	; TODO: Perform some checks!

	; Assign the values
	for i = 0, n_elements(id)-1 do begin
	    indices[i].id = id[i]
	    indices[i].name = strcompress(name[i], /remove_all)
	    indices[i].units = strcompress(units[i], /remove_all)
	    indices[i].passband=[pass0[i],pass1[i]]
	    indices[i].blue_cont=[blcnt0[i],blcnt1[i]]
	    indices[i].red_cont=[rdcnt0[i],rdcnt1[i]]
	endfor

	return, indices
END

