;+
; NAME:
;       MDAP_NOISE_FROM_RESIDUAL
;
; PURPOSE:
;       Calculate the error vector based on the residual of a fitted vector with
;       respect to the data.
;
; CALLING SEQUENCE:
;       MDAP_NOISE_FROM_RESIDUAL, data, mask, model, sigma, bin_width=bin_width, $
;                                 min_width=min_width, /roaming
;
; INPUTS:
;       data dblarr[D]
;               Data values for D points.
;       
;       mask dblarr[D]
;               Bad pixel mask for the data: 0-good, 1-bad.  Changed upon output
;               to reflect minimum number of unmasked pixels necessary for sigma
;               calculation
;       
;       model dblarr[D]
;               Model values
;
; OPTIONAL INPUTS:
;       bin_width integer
;               Number of pixels to include in estimation of noise at any given
;               position in the dataset.
;
;               IF NOT PROVIDED, the default it to assume a constant error for
;               all data by using all data to calculate the rms noise.
;
;       min_width integer
;               Minimum number of pixels required for determining the standard
;               deviation in a bin.
;
;               IF NOT PROVIDED, the default is (arbitrarily) set to 5 pixels. 
;
; OPTIONAL KEYWORDS:
;       /roaming
;               If set, the bin used to calculate the rms noise is stepped by a
;               single pixel instead of the full bin width.
;
; OUTPUT:
;       sigma dblarr[D]
;               Error in the data determined by the residual of the data with
;               respect to the model.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Define how the edges should be treated
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       03 Nov 2014: (KBW) Based on the algorithm in MDAP_MEASURE_INDICES by L.
;                    Coccato (PDAP v0_8)
;       04 Nov 2014: (KBW) Impose a minimum bin width
;-
;------------------------------------------------------------------------------

PRO MDAP_NOISE_FROM_RESIDUAL, $
                data, mask, model, sigma, bin_width=bin_width, min_width=min_width, roaming=roaming

        ndata = n_elements(data)                        ; Number of data points
        diff = data-model                               ; Get the residual vector

        ; There is only one bin (one calculation of sigma) if ...
        one_bin = n_elements(bin_width) gt 0 ? 0 : 1    ; bin_width is not provided
;       if one_bin eq 0 and bin_width gt ndata then $   ; or it's larger than # of data points
        if one_bin eq 0 && bin_width gt ndata then $   ; or it's larger than # of data points
            one_bin = 1

        if one_bin eq 1 then begin                      ; Only one sigma for all data values
            diff = data-model                           ; Get the residual vector
            indx = where(mask lt 0.5)                   ; Find the unmasked pixels

            ; Calculate sigma and define the array
            sigma = make_array(ndata, /double, value=robust_sigma(diff[indx]))
            return                                      ; DONE
        endif

        if n_elements(min_width) eq 0 then $            ; Set the default minimum width
            min_width = 5

        if bin_width lt min_width then $                ; Check that the input bin width is valid
            message, 'Minimum bin width is '+MDAP_STC(min_width,/integer)+'!'

        if keyword_set(roaming) then begin      ; Bin roams 1 pixel at a time
            nbins = ndata                       ; Number of bins same as the number of data points
            bin_step = 1                        ; Bin step is 1 pixel
            first_bin = 0                       ; First bin center is at the first pixel
        endif else begin
            nbins = fix(ndata/bin_width)        ; Bins are independent
            bin_step=bin_width                  ; Bin step is the same as the bin width
            first_bin=bin_width/2               ; First bin is put at the edge of the data
        endelse

        data_center = findgen(ndata)            ; Define the data indices
        sigma_ = dblarr(nbins)                  ; Array to hold the sigma values
        bin_center = dblarr(nbins)              ; ... and the bin centers
        mask_ = dblarr(nbins)                   ; ... and the updated mask

        for i=0,nbins-1 do begin
            bin_center[i] = first_bin+i*bin_step        ; Set the bin center

            ; Pixels to include are those within the bin and not masked.  When
            ; "roaming," the first and last bin_width/2 pixels will not be based
            ; on a full bin.
            indx = where( data_center ge bin_center[i]-bin_width/2 and $
                          data_center lt bin_center+bin_width/2 and mask lt 0.5 )

            if n_elements(indx) lt min_width then begin ; Require at least min_width pixels!
                mask_[i] = 1.0d                         ; Mask the pixel
            endif else $
                sigma_[i] = robust_sigma(diff[indx])    ; Calculate sigma
        endfor

        if keyword_set(roaming) then begin
            sigma = sigma_                      ; Sigma vector generated via the above identically
            mask = mask_                        ; ... as is the mask
            indx = where(mask gt 0.5)
            if indx[0] ne -1 then $
                sigma[indx] = 1.0d              ; Force nominal error for masked pixels
            return
        endif

        ; Replace masked pixels with nearest unmasked pixel to minimize effects
        ; in interpolation
        indx = where(mask_ gt 0.5, complement=nindx)
        if nindx[0] eq -1 then $
            message, 'All pixels masked!'
        if indx[0] ne -1 then begin
            nmask = n_elements(indx)
            for i=0,nmask-1 do begin
                nindx_ = (sort(abs(nindx-indx[i])))[0]
                sigma_[indx[i]] = sigma_[nindx[nindx_]]
            endfor
        endif

        sigma = interpol(sigma_, bin_center, data_center)       ; Interpolate the sigma
        mask = interpol(mask_, bin_center, data_center)         ; Interpolate the mask

        indx = where(mask lt 0.5, complement=nindx)             ; Force binary mask
        mask[indx] = 0.0d
        mask[nindx] = 1.0d
        sigma[nindx] = 1.0d                                     ; Force nominal error

END

;    nbins = fix(n_elements(spc)/20)
;    bins = mdap_range(0, n_elements(spc),nbins+1)
;    fake_x = findgen(n_elements(spc))
;    rms_per_bin = fltarr(nbins)
;    bin_center = fltarr(nbins)
;     
;    for kk = 0, nbins-1 do begin
;       indici = where(fake_x ge bins[kk] and fake_x lt bins[kk+1])
;       rms_per_bin[kk] = stddev(res[indici])  
;       bin_center[kk] = (bins[kk+1]+ bins[kk])*0.5
;    endfor
;    noise = interpol(rms_per_bin,bin_center,findgen(n_elements(spc)))

