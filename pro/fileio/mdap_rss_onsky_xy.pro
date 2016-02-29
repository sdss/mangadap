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
;       24 Jul 2015: (KBW) DRP convention is that x offsets are positive
;                          toward the west and negative toward the east,
;                          which is opposite to what is currently done
;                          in MDAP_CUBE_ONSKY_XY for the CUBE files.
;                          I've switched this convension to be
;                          consistent with the radial binning approach.
;       26 Oct 2015: (KBW) Apply offset to coordinates based on
;                          difference between IFURA/IFUDEC and
;                          OBJRA/OBJDEC.  AS OF DRP VERSION 1_5_1, the
;                          convention is that IFURA/IFUDEC is the
;                          nominal center/pointing of the IFU dither
;                          patter, whereas OBJRA/OBJDEC is the nominal
;                          center of the galaxy.  The two coordinates
;                          are most often the same, except for
;                          intentional offsets (like 7443-12701).
;       29 Feb 2016: (KBW) Fixed cos() argument to be in radians!
;-
;-----------------------------------------------------------------------

PRO MDAP_RSS_ONSKY_XY, $
        file, header, skyx, skyy

        MDAP_READ_FITS_EXTENSION, file, 'XPOS', skyx, /to_double
        MDAP_READ_FITS_EXTENSION, file, 'YPOS', skyy, /to_double

        skyx=transpose(temporary(-skyx))            ; Transpose the matrices to match spectra
        skyy=transpose(temporary(skyy))

        objra = SXPAR(header, 'OBJRA')              ; Object RA in degrees
        objdec = SXPAR(header, 'OBJDEC')            ; Object DEC in degrees
        
        ifura = SXPAR(header, 'IFURA')              ; IFU RA in degrees
        ifudec = SXPAR(header, 'IFUDEC')            ; IFU DEC in degrees

        x_offset = (ifura-objra)*cos(objdec*!PI/180.)*3600. ; Offset in RA
        y_offset = (ifudec-objdec)*3600.           ; Offset in DEC

        print, 'IFU RA/DEC:', ifura, ifudec
        print, 'OBJ RA/DEC:', objra, objdec
        print, 'Applying offset (arcsec): ', x_offset, y_offset

        skyx = skyx + x_offset
        skyy = skyy + y_offset

END


