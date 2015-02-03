;+
; NAME:
;       MDAP_RSS_ONSKY_XY
;
; PURPOSE:
;       Set the on-sky x and y coordinates of the RSS files using data from the
;       fits file.  The order of the arrays are rearranged to match the
;       restructured specta (see mdap_restructure_rss.pro).
;
; CALLING SEQUENCE:
;       MDAP_RSS_ONSKY_XY, file, skyx, skyy
;
; INPUTS:
;       file string
;               Fits file with x and y coordinates in extensions 7 and 8.
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
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       09 Sep 2014: (KBW) Original Implementation
;       01 Feb 2015: (KBW) Switch to reading the file using the name of
;                          the extension instead of its number
;-
;-----------------------------------------------------------------------

PRO MDAP_RSS_ONSKY_XY, $
        file, skyx, skyy

        MDAP_READ_FITS_EXTENSION, file, 'XPOS', skyx, /to_double
        MDAP_READ_FITS_EXTENSION, file, 'YPOS', skyy, /to_double

        skyx=transpose(temporary(skyx))                 ; Transpose the matrices to match spectra
        skyy=transpose(temporary(skyy))
END
        


