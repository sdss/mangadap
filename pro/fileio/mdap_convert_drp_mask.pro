;+
; NAME:
;       MDAP_CONVERT_DRP_MASK
;
; PURPOSE:
;       Convert the mask read from the DRP MASK extension to a binary
;       unmasked/masked array for use in the DAP routines.
;
; CALLING SEQUENCE:
;       MDAP_CONVERT_DRP_MASK, mask
;
; INPUTS:
;       mask dblarr[N][T]
;               On input contains the maskbits read from the MASK
;               extension of a DRP output file.  On output, is a
;               bad-pixel mask for the DRP FLUX extension, with 0.
;               meaning unmasked, 1. meaning masked.
;
;               WARNING: If no maskbits exist, the return array will
;               consist of all unmasked pixels!
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       !! mask is replaced upon output !!
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
;       MDAP_DEFINE_BAD_PIXEL_FLAGS
;
; REVISION HISTORY:
;       06 Feb 2015: Original implementation by K. Westfall
;-
;-----------------------------------------------------------------------

PRO MDAP_DEFINE_BAD_PIXEL_FLAGS, $
                maskbit_type, flag_name, exists

        ; Define the maskbit type
        maskbit_type = 'MANGA_DRPPIXFLAG'

        ; Define ALL the maskbits to set to a bad pixel
        flag_name = [ 'LOWCOV', 'NOCOV' ]
;       flag_name = [ 'DEADFIBER', 'BADDATA', 'LOWCOV', 'NOCOV' ]

        ; Determine which maskbits are actually defined, and warn the
        ; user if any of the above are not defined
        nn = n_elements(flag_name)
        exists = intarr(nn)
        for i=0,nn-1 do begin
            exists[i] = SDSS_FLAGEXIST(maskbit_type, flag_name[i])
            if exists[i] ne 1 then $
                print, 'WARNING: No '+maskbit_type+' maskbit='+flagname[i]+' has been defined!'
        endfor
END


PRO MDAP_CONVERT_DRP_MASK, $
                mask

        mask_ = mask                                            ; Save the input

        sz = size(mask)                                         ; Get the input dimensions
        dim = intarr(sz[0])
        for i=0,sz[0]-1 do $
            dim[i] = sz[i+1]

        ; Reset the mask to 0.0 while keeping the same input dimensions
        mask = make_array(dimension=dim, /byte, value=0)

        ; Get the maskbits to flag as bad pixels
        MDAP_DEFINE_BAD_PIXEL_FLAGS, maskbit_type, flag_name, exists

        ; If none of the maskbits exist, warn the user and return;
        ; output mask will be all 0.
        indx = where(exists gt 0, count)
        if count eq 0 then begin
            print, 'WARNING: No maskbits defined!'
            mask = double(mask)                         ; Convert to double before returning
            return
        endif

        ; Mask all the pixels defined to be bad
        for i=0,count-1 do $
            mask = mask or ((mask_ and SDSS_FLAGVAL(maskbit_type, flag_name[indx[i]])) ne 0)

        ; Convert to double to use with the rest of the DAP
        mask = double(mask)

;       bad_indx_drp = where(mask_ gt 0)
;       bad_indx_dap = where(mask gt 0)

;       print, n_elements(mask), n_elements(bad_indx_drp), n_elements(bad_indx_dap)
END

