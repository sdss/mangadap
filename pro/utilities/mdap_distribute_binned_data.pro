;+
; NAME:
;       MDAP_DISTRIBUTE_BINNED_DATA
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_DISTRIBUTE_BINNED_DATA, binned_data, binned_index, distributed_data, llim=llim, $
;                                    ulim=ulim, default=default
;
; INPUTS:
;       binned_data dblarr[B]
;               Data associated with a binned spectrum.
;
;       binned_indx intarr[N]
;               Array giving the index of the bin in which each original
;               spectrum was placed.
;
; OPTIONAL INPUTS:
;       llim double
;               Lower limit allowed for output value
;
;       ulim double
;               Upper limit allowed for output value
;
;       default double
;               Default value to use if spectrum was not part of any bin
;               or if the value falls outside of the allowed range.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       distributed_data dblarr[N]
;               A redistribution of the binned data to those spectra in
;               the bin.
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
;       04 Dec 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_DISTRIBUTE_BINNED_DATA, $
                binned_data, binned_index, distributed_data, llim=llim, ulim=ulim, default=default

        ns = n_elements(binned_index)       ; Number of spectra
        nb = n_elements(binned_data)        ; Number of bins

        if n_elements(default) eq 0 then $
            default = 0.0d

        distributed_data = make_array(ns, /double, value=default)

        for i=0,ns-1 do begin
            if binned_index lt 0 then $     ; Spectrum not part of a bin
                continue
            distributed_data[i] = binned_data[binned_index[i]]
        endfor

        if n_elements(llim) ne 0 then begin
            indx = where(distributed_data[i] lt llim)
            if indx[0] ne -1 then
                distributed_data[indx] = default
        endif

        if n_elements(ulim) ne 0 then begin
            indx = where(distributed_data[i] gt ulim)
            if indx[0] ne -1 then
                distributed_data[indx] = default
        endif
END


