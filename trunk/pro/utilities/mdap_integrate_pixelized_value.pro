;+
; NAME:
;       MDAP_INTEGRATE_PIXELIZED_VALUE
;
; PURPOSE:
;       Performs the integral:
;
;               integral = \int_xrange[0]^xrange[1] y dx
;
;       using a Riemann sum, where x is expected to be at the center of
;       the pixel.  The error is estimated by a simple propagation of
;       the error (i.e. error in the sum).
;
; CALLING SEQUENCE:
;       MDAP_INTEGRATE_PIXELIZED_VALUE, x, y, ye, mask, xrange, integral, integral_err, err=err, $
;                                       /geometric
;
; INPUTS:
;       x dblarr[]
;               X coordinate vector.  It is expected to be sorted and
;               contiguous; x coordinates are assumed to be linearly or
;               geometrically spaced based on the input keyword
;               (/geometric); the provided coordinate is expected to be
;               at the (geometric) center of the pixel.
;
;       y dblarr[]
;       ye dblarr[]
;               Y coordinate vector and its error.  Samples of the
;               function to integrate at the provided x coordinates.
;
;       mask dblarr[]
;               Bad pixel mask (0-good;1-bad) for the input function.
;               Masked pixels are ignored in the integral.  If not
;               pixels in the provided xrange are good, the integral and
;               its error are both returned as 0.  However, the error
;               value is set to 0!
;
;       xrange dblarr[2]
;               The start (index 0) and end (index 1) of the range over
;               which to integrate y.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /geometric
;               The sampling of the wavelength vector is geometric
;               (e.g., logarithmic); required by MDAP_BIN_EDGES to set
;               the edges of the wavelength pixels.
;
; OUTPUT:
;       integral double
;       integral_err double
;               The integral of function and its error.
;
; OPTIONAL OUTPUT:
;       err integer
;               Flag that an error occurred during the calculation.  0
;               for no error, 1 for error.  Currently only returned as
;               err=1 if the xrange does not overlap the provided x
;               vector.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Errors seem bogus, but following nominal error propagation
;         procedure...
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_BIN_EDGES()
;
; REVISION HISTORY:
;       17 Aug 2015: Moved from MDAP_GET_SPECTRAL_INDEX by K. Westfall (KBW), 
;                    for use in other routines, e.g.
;                    MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS
;       26 Oct 2015: (KBW) Check that the spectrum to integrate over has
;                          pixels that are non-zero.  Otherwise, return
;                          err=1.
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Bin edge [0] is the lower edge of pixel x[0], the upper edge of pixel x[0] is
; bin edge [1], etc.  The values of x are expected to be linearly or
; geometrically spaced.
FUNCTION MDAP_BIN_EDGES, $
                x, geometric=geometric

        ; Get the number of bin centers
        n = n_elements(x)
        if n eq 1 then $
            message, 'Cannot determine edges of a single bin!'

        ; Convert the bin centers to their logarithm if the bins are
        ; logarithmically spaced
        xx = x
        if keyword_set(geometric) then $
            xx = alog(x)

        ; Calculate the n+1 edges of the bins
        edges = ([ xx[0]-(xx[1]-xx[0])/2.0d , xx[0:n-2] + (xx[1:n-1]-xx[0:n-2])/2.0, $
                   xx[n-1]+(xx[n-1]-xx[n-2])/2.0d ])

        ; Convert them back to linear space, if necessary
        if keyword_set(geometric) then $
            edges = exp(edges)

        ; Return the bin edges
        return, edges
END

;-------------------------------------------------------------------------------
; Performs the integral:
;
;       integral = \int_xrange[0]^xrange[1] y dx
;
; using a Riemann sum, where x is expected to be at the center of the pixel.
; The error is estimated by a simple propagation of the error (i.e. error in the
; sum).
;
; x is expected to be sorted and contiguous; x coordinates are assumed to be
; linearly or geometrically spaced; the provided coordinate is expected to be at
; the (geometric) center of the pixel.
PRO MDAP_INTEGRATE_PIXELIZED_VALUE, $
                x, y, ye, mask, xrange, integral, integral_err, err=err, geometric=geometric

        err=0                                                   ; Initialize error
        n = n_elements(x)                                       ; Number of pixels
        bin_edges = MDAP_BIN_EDGES(x, geometric=geometric)      ; The n+1 edges of the n pixels
        if xrange[1] lt bin_edges[0] || xrange[0] gt bin_edges[n] then begin    ; No overlap
;           print, xrange[0], xrange[1], bin_edges[0], bin_edges[n]
            integral = 0.0d
            integral_err = 1.0d
            err=1                                               ; Return an error
            return
        endif

        ; Set the values to 0.0 for masked pixels
        indx = where(mask gt 0.5, count)
        ym = y
        yme = ye
        if count ne 0 then begin
            ym[indx] = 0.0d
            yme[indx] = 0.0d
        endif

        ; Interval is smaller than the nearest pixel
        indx = where(bin_edges gt xrange[0] and bin_edges lt xrange[1], ni)
        if ni eq 0 then begin
            indx = where(bin_edges gt xrange[1])
            integral = ym[indx[0]-1]*(xrange[1]-xrange[0])
            integral_err = yme[indx[0]-1]*(xrange[1]-xrange[0])
;            print, 'Interval smaller than a pixel!'
;            print, integral, integral_err
            return
        endif

        ; Initialize output
        integral = 0.0d
        integral_err = 0.0d

        ; Check that there are non-zero values to integrate
        nonzero = where(ym[indx] ne 0, nnz)
        if nnz eq 0 then begin
;            print, 'No non-zero fluxes!'
            integral = 0.0d
            integral_err = 1.0d
            err=1                                               ; Return an error
            return
        endif

        ; Add partial pixel at the blue end
        if indx[0] gt 0 then begin
            integral = integral + ym[indx[0]-1]*(bin_edges[indx[0]]-xrange[0])
            integral_err = integral_err + (yme[indx[0]-1]*(bin_edges[indx[0]]-xrange[0]))^2
        endif

        ; Add full pixels
        if ni ge 2 then begin
            integral = integral + total( ym[indx[0:ni-2]]*(bin_edges[indx[1:ni-1]] $
                                                          - bin_edges[indx[0:ni-2]]) )
            integral_err = integral_err + total( (yme[indx[0:ni-2]]*(bin_edges[indx[1:ni-1]] $
                                                                     - bin_edges[indx[0:ni-2]]))^2 )
        endif

        ; Add partial pixel at the red end
        if indx[ni-1] lt n then begin
            integral = integral + ym[indx[ni-1]]*(xrange[1]-bin_edges[indx[ni-1]])
            integral_err = integral_err + ( yme[indx[ni-1]]*(xrange[1]-bin_edges[indx[ni-1]]) )^2
        endif

        integral_err = sqrt(integral_err)
;       if integral eq 0 || integral_err eq 0 then begin
;           print, ni, integral, integral_err
;           plot, x, ym
;           oplot, x[indx], ym[indx], color=200
;           print, x[indx]
;           print, ym[indx]
;           print, y[indx]
;           stop
;       endif
END

