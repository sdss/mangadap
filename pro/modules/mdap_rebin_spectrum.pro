;+
; NAME:
;   rebin_spectrum
;
; PURPOSE:
;   Rebin a 1-D spectrum using CIC assignment onto arbitrary pixel boundaries.
;
; CALLING SEQUENCE:
;   yflux = rebin_spectrum( xflux, xwave, ywave )
;
; INPUTS:
;   xflux      - Flux vector [NOLD]
;   xwave      - Wavelengths of pixel edges for input XFLUX [NOLD+1],
;                or the central wavelength of each pixel [NOLD]
;   ywave      - Wavelengths of pixel edges for output YFLUX [NNEW+1],
;                or the central wavelength of each pixel [NNEW]
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   yflux      - Rebinned spectrum
;
; COMMENTS:
;   This function does a straight cloud-in-cell re-assignment of flux
;   from one spectrum to another.  The boundaries of the flux in pixel #i,
;   XFLUX[i], in the first spectrum is assumed to be uniformly distributed
;   between the wavelengths [XWAVE[i],XWAVE[i+1]].  This flux is re-assigned
;   to YFLUX, whose pixel #j is assumed to cover [YWAVE[j],YWAVE[j+1]].
;   This algorithm is strictly flux-conserving for the wavelengths that
;   overlap.
;
;   If one wavelength grid is an integral multiple of the other, than this
;   is equivalent to using the IDL REBIN command.  For example, the following
;   exactly puts the flux from the XFLUX spectrum into bins twice as big:
;     IDL> xflux=randomu(1234,100)
;     IDL> xwave=findgen(101)
;     IDL> ywave=findgen(51)*2
;     IDL> yflux=rebin_spectrum(xflux,xwave,ywave)
;   This could also be accomplished with
;     IDL> yflux2=rebin(xflux,50)
;   In this example, the two pixels of the input spectrum span the
;   wavelengths [0,1] and [1,2].  The first pixel of the output spectrum
;   spans [0,2].
;
;   Both XWAVE and YWAVE must be in strictly ascending order, and cannot
;   repeat any wavelengths within those vectors.
;
;   One can also call this routine with XWAVE having the same number of
;   elements as XFLUX.  In that case, we assume that the wavelengths are
;   at the center of each pixel, and internal to this routine compute
;   (by linear interpolation) where the pixel boundaries are.  For this case,
;   we also interpret YWAVE as being at the pixel centers.
;
; BUGS:
;   This should probably be implmented in C for better speed.
;   We do not test that XWAVE and YWAVE are strictly ascending.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   18-Jan-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function mdap_rebin_spectrum, xflux, xwave, ywave

   if (n_params() LT 3) then begin
      print, 'Syntax - yflux = rebin_spectrum(xflux, xwave, ywave)'
      return, 0
   endif
   npts = n_elements(xflux)
   if (n_elements(xwave) EQ npts) then begin
      ; In the case where the central wavelengths are specified,
      ; we need to compute the wavelengths at the edges of the pixels
      xtmp = interpol(xwave, 2*lindgen(npts)+1, 2*lindgen(npts+1))
      nnew = n_elements(ywave)
      ytmp = interpol(ywave, 2*lindgen(nnew)+1, 2*lindgen(nnew+1))
      return, mdap_rebin_spectrum(xflux, xtmp, ytmp)
   endif
   if (n_elements(xwave) NE npts+1) then $
    message, 'Number of elements in XWAVE must equal N or N+1'

   nnew = n_elements(ywave)-1
   yflux = make_array(nnew, type=size(xflux,/type)) ; Same data type as input

   ;----------
   ; If no overlapping wavelength region, then exit without assigning
   ; any flux.

   if (xwave[npts] LT ywave[0] OR xwave[0] GT ywave[nnew]) then return, yflux

   ;----------
   ; Start by skipping any pixels in the input spectrum where there is
   ; no place to assign the flux

   j1 = 0L
   i = 0L
   while (xwave[i] LT ywave[0]) do i = i + 1
   ; Assign some fraction of the flux from the first pixel...
   if (i GT 0) then begin
      dflux = xflux[i-1] / (xwave[i] - xwave[i-1])
      while (ywave[j1+1] LT xwave[i]) do begin
         yflux[j1] += dflux * (ywave[j1+1] - ywave[j1])
         j1 = j1 + 1
      endwhile
      yflux[j1] += dflux * ((xwave[i] < ywave[j1+1]) - ywave[j1])
   endif
   j2 = j1

   ;----------
   ; Loop over each pixel in the input spectrum, and assign its flux

   while (i LT npts) do begin
      while (ywave[j1+1] LT xwave[i] AND j1 LT nnew-1) do j1 = j1+1
      while (ywave[j2+1] LT xwave[i+1] AND j2 LT nnew-1) do j2 = j2+1

      dflux = xflux[i] / (xwave[i+1] - xwave[i])
      if (j1 EQ j2) then begin
         yflux[j1] += dflux * ((ywave[j1+1] < xwave[i+1]) - xwave[i])
      endif else begin
         yflux[j1] += dflux * (ywave[j1+1] - xwave[i])
         for k=j1+1, j2-1 do $
          yflux[k] += dflux * (ywave[k+1] - ywave[k])
         if (ywave[j2+1] GE xwave[i+1]) then $
          yflux[j2] += dflux * (xwave[i+1] - ywave[j2])
      endelse
      i = i + 1
      if (ywave[j2+1] LT xwave[i]) then i = npts
   endwhile

   return, yflux
end
;------------------------------------------------------------------------------
