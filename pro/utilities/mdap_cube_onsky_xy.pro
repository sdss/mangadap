;+
; NAME:
;       MDAP_CUBE_ONSKY_XY
;
; PURPOSE:
;       Set the on-sky x and y coordinates of the data cube using the header WCS
;       solution.
;
; CALLING SEQUENCE:
;       MDAP_CUBE_ONSKY_XY, header, skyx, skyy
;
; INPUTS:
;       header hdu
;               Header of the DRP-produced fits file, obtained from HEADFITS().
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       skyx dblarr[][]
;               NAXIS1*NAXIS2 by T array with the x-coordinates of each pixel,
;               where T is the number of spectral channels.
;
;       skyy dblarr[]][]
;               NAXIS1*NAXIS2 by T array with the y-coordinates of each pixel,
;               where T is the number of spectral channels.
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
;       SXPAR()
;       XYAD
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       09 Sep 2014: (KBW) Original Implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_CUBE_ONSKY_XY, $
        header, skyx, skyy

        ; TODO: Does not check if these header keywords exist!
        nx = SXPAR(header, 'NAXIS1')                    ; Number of x-coordinates
        ny = SXPAR(header, 'NAXIS2')                    ; Number of y-coordinates
        ns = SXPAR(header, 'NAXIS3')                    ; Number of spectral channels

        x = MAKE_ARRAY(nx*ny, /double)                  ; Allocate the arrays
        y = MAKE_ARRAY(nx*ny, /double)
        for i=0,nx-1 do begin
            for j=0,ny-1 do begin
                x[i*ny+j,*] = i+1                       ; Set the pixel coordinates
                y[i*ny+j,*] = j+1
            endfor
        endfor

        XYAD, header, x, y, skyx_, skyy_                ; Get the coordinates

        skyx = make_array(nx*ny, ns, /double)
        skyy = make_array(nx*ny, ns, /double)
        for i=0,ns-1 do begin
            skyx[*,i] = skyx_
            skyy[*,i] = skyy_
        endfor
END
        


