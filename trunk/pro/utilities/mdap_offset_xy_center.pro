;+
; NAME:
;       MDAP_OFFSET_XY_CENTER
;
; PURPOSE:
;       Compute the distance of each spaxel from the center of the targetted
;       object using the keywords 'OBJRA' and 'OBJDEC' from the supplied header.
;       As this is a calculation of distance, the fomula includes the
;       approximate correction for the cosine of the declination.
;
; CALLING SEQUENCE:
;       MDAP_OFFSET_XY_CENTER, header, skyx, skyy
;
; INPUTS:
;       header hdu
;               Header read from the DRP-produced fits file (CUBE or RSS).
;
;       skyx dblarr[]
;               On-sky x coordinates of the spectra.  This is converted to x
;               distance from the center.
;
;       skyy dblarr[]
;               On-sky y coordinates of the spectra.  This is converted to y
;               distance from the center.
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
;       - Check that other procedures expect input that has the distance from
;         the center!
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       09 Sep 2014: (KBW) Original Implementation
;       15 Sep 2014: (KBW) Just offsets to the center using the WCS units
;-
;------------------------------------------------------------------------------

PRO MDAP_OFFSET_XY_CENTER, $
                header, skyx, skyy

        ; TODO: Does not check that keywords exist!
        targx = SXPAR(header, 'OBJRA')                  ; Read the target RA ...
        targy = SXPAR(header, 'OBJDEC')                 ; ... and declination

        ; TODO: Eventually pull this from DRPall; TARGETRA, TARGETDEC
        ;       OBJRA is the same as IFURA, etc

        skyx = (skyx-targx)*cos(targy*!DtoR)    ; Distance with respect to the target center
        skyy = (skyy-targy)

END

