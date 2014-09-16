;+
; NAME:
;	MDAP_RESTRUCTURE_CUBE
;
; PURPOSE:
;	Restructure the input data cube into a 2D array for use with main DAP
;	procedures.  The input data arrays are expected to be 3D, [N][M][T],
;	where N is along the on-sky X axis, M is along the on-sky Y axis, and T
;	is along the spectral dimension.  This cube is put into a 2D array,
;	[N*M][T], with the array ordered such that the (n,m)th spectrum is in
;	index [n*NY+m,*] of the array.
;	
; CALLING SEQUENCE:
;	MDAP_RESTRUCTURE_CUBE, flux, ivar, mask
;
; INPUTS:
;	flux dblarr[N][M][T]
;		Input cube read from a DRP-produced 'CUBE' file.  Converted to
;		new format on output.
;
;	ivar dblarr[N][M][T]
;		Input inverse-variance cube read from a DRP-produced 'CUBE'
;		file.  Converted to new format on output.
;
;	mask dblarr[N][M][T]
;		Input pixel mask cube read from a DRP-produced 'CUBE' file.
;		Converted to new format on output.
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
; TODO:
;	- THIS NEEDS TO BE FASTER!
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	09 Sep 2014: (KBW) Original implementation
;	09 Sep 2014: (KBW) Include mask data
;-
;------------------------------------------------------------------------------

PRO MDAP_RESTRUCTURE_CUBE, $
		flux, ivar, mask

	sz=size(flux)
	nx = sz[1]					; Get the size in the x-dimension
	ny = sz[2]					; Get the size in the y-dimension
	ns = sz[3]					; Number of spectral channels

	flux_ = flux					; Copy the input
	ivar_ = ivar
	mask_ = mask

	flux=dblarr(nx*ny, ns, /nozero)			; Initialize the arrays (quickly)
	ivar=dblarr(nx*ny, ns, /nozero)
	mask=dblarr(nx*ny, ns, /nozero)
	for i=0,nx-1 do begin
	    for j=0,ny-1 do begin
		flux[i*ny+j,*] = flux_[i,j,*]		; Reorganize the data
		ivar[i*ny+j,*] = ivar_[i,j,*]
		mask[i*ny+j,*] = mask_[i,j,*]
	    endfor
	endfor

END



