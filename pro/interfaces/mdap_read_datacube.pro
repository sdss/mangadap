;+
; NAME:
;   MDAP_READ_DATACUBE
;
; PURPOSE:
;	Convert RSS and datacube output from DRP to vectors that are processed
;	by other DAP modules.
;
; CALLING SEQUENCE:
;	MDAP_READ_DATACUBE,$
;		datacube_name, data, error, wavelength, x2d, y2d, signal, noise, cdelt1, cdelt2,$
;		header2d, lrange=lrange, keep_original_step=keep_original_step,$
;		x2d_reconstructed=x2d_reconstructed, y2d_reconstructed=y2d_reconstructed,$
;		signal2d_reconstructed=signal2d_reconstructed,$
;		noise2d_reconstructed=noise2d_reconstructed, version=version,$
;		number_of_fibres=number_of_fibres, use_total=use_total
;
; INPUTS:
;   datacube_name [string]
;
;		Name of the fits file containing the input in datacube or RSS
;		formats.
;		
;		For inputs in datacube format, the file must be a fits
;		file with the following extensions:
;
;		    1. flux in 1e-17 ergs/s/cm^{-2}/angstrom. This extension
;		    must be a 3D array, with the wavelength direction along the
;		    3rd axis.
;
;		    2. Inverse variance associated to the first extension
;
;		    3. Wavelength solution (1D dbl array). Constant linear step
;		    (preferred)? Constant Log-step? Constant ln-step?
;		    
;		    4. RSS format (row-stacked spectra) (NOT USED) 
;		    
;		    5. coverage map; (NOT USED)
; 
;		For inputs in RSS format, the file must be a fits file with the
;		following extensions:
;
;		    1. flux in 1e-17 ergs/s/cm^{-2}/angstrom. This extension
;		    must be a 3D array, with the wavelength direction along the
;		    3rd axis.
;		    
;		    2. Inverse variance associated to the first extension
;		    
;		    3. Wavelength solution (1D dbl array). Constant linear step
;		    (preferred)? Constant Log-step? Constant ln-step?
;		    
;		    4. X position table.  Since position is a function of
;		    wavelength this is an array of the same size as the flux
;		    array.  Coordinates should be in arcseconds relative to the
;		    center of the data cube (i.e., roughly the center of the
;		    galaxy).
;		    
;		    5. Y position table.
;
; OPTIONAL INPUTS:
;    lrange vector[2]
;		Wavelength range (in angstrom) used in signtal-to-noise calculation.
;		Default: use the entire spectra range.
;
; OPTIONAL KEYWORDS:
;    /keep_original_step
;		If set, the wavelength output vector will be the same as the one
;		define from the input fits file. The default is to re-construct
;		(and  therefore re-inrpolate the galaxy and error spectra) the
;		output wavelength vector with constant ang/pixel step using the
;		minimum ang/pixel step that is stored in the wavelength
;		solution. For MANGA, it is suggested to set this keyword.
;
;    /use_total
;		If set, the signal is the total of the counts in the selected
;		wavelength range, the noise is the sum in quadrature of the
;		noises in the selected range. Useful for emission lines.
;
; OUTPUT:
;    data dblarr[N][M][T] or dblarr[N][T]
;		Datacube of galaxy spectra. Spectra are resampled over the
;		vector wavelength.  N and M are spatial dimensions, T is the
;		spectral dimension.  RSS spectra only have one spatial
;		dimension.
;
;    error dblarr[N][M][T] or dblarr[N][T]
;		Datacube of error in galaxy spectra, returned by DRP. N and M
;		are spatial dimensions, T is the spectral dimension.  RSS
;		spectra only have one spatial dimension.
;
;    wavelength dblarr[T]
;		Wavelength coordinates of each pixel. The vector is constructed
;		with constant linear step in ang/pixel (unless the
;		/keep_original_step keyword is selected).  If input spectra have
;		a logarithmic sampling, the minimum available step is used (e.g.
;		log_lam[1]-log_lam[0]).
;
;    x2d dblarr[N][M] or dblarr[N]
;		Spatial X coordinate in arcseconds; 0 is the center of the field
;		of view.
;
;    y2d dblarr[N][M] or dblarr[N]
;		Spatial Y coordinate in arcseconds; 0 is the center of the field
;		of view.
;
;    signal dblarr[N][M] or dblarr[N]
;		Mean galaxy signal per angstrom, obtained considering the full
;		wavelength range (or only the range specified by lrange).
;		Calculation uses original spectra, not those resampled over the
;		vector wavelength.  Signal calculation uses median(signal),
;		unless the keyword /use_total is set.
;
;    noise dblarr[N][M] or dblarr[N]
;		Mean galaxy noise per angstrom, obtained considering the full
;		wavelength range (or only the range specified by lrange).
;		Calculation uses original spectra, not those resampled over the
;		vector wavelength.  Noise calculation uses median(noise), unless
;		the keyword /use_total is set.
;
;    cdelt1 double
;		Spatial sampling along x direction (arcsec/pixel).  Set to 0.5
;		for RSS input.
;
;    cdelt2 double
;		Spatial sampling along y direction (arcsec/pixel).  Set to 0.5
;		for RSS input.
;
;    header2d   string[]
;		The header for two-dimensional fields.  TODO: RSS format output?
;
; OPTIONAL OUTPUT:
;    x2d_reconstructed dblarr[N][M]
;		For RSS input, the x2d coordinates are resampled over a 2D map
;		with fixed 0.5 arcsec/pixel sampling.  For datacube input, this
;		is exactly x2d.
;
;    y2d_reconstructed dblarr[N][M]
;		For RSS input, the y2d coordinates are resampled over a 2D map
;		with fixed 0.5 arcsec/pixel sampling.  For datacube input, this
;		is exactly y2d.
;
;    signal2d_reconstructed dblarr[N][M]
;		For RSS input, the signal values are resampled over a 2D map
;		with fixed 0.5 arcsec/pixel sampling.  For datacube input, this
;		is exactly signal.
;
;    noise2d_reconstructed dblarr[N][M]
;		For RSS input, the noise values are resampled over a 2D map
;		with fixed 0.5 arcsec/pixel sampling.  For datacube input, this
;		is exactly noise.
;
;    version string
;		Module version. If requested, code simply returns the version
;		flag.
;
;
; TODO:
;	- header output for RSS data?
;	- possible to force only ONE output size for both datacube and RSS data?
;	- split into two programs; one to read RSS another to read CUBE, both of
;	  which produce the same final output
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
;	01 Sep 2014: Copied from v0_8 version by L. Coccato.
;	01 Sep 2014: (KBW) Formating and edits to deal with DRP v1_0_0 header
;	information
;	05 Sep 2014: (KBW) More formatting
;-
;------------------------------------------------------------------------------

PRO MDAP_READ_DATACUBE,$
		datacube_name, data, error, wavelength, x2d, y2d, signal, noise, cdelt1, cdelt2,$
		header2d, lrange=lrange, keep_original_step=keep_original_step,$
		x2d_reconstructed=x2d_reconstructed, y2d_reconstructed=y2d_reconstructed,$
		signal2d_reconstructed=signal2d_reconstructed,$
		noise2d_reconstructed=noise2d_reconstructed, version=version,$
		number_of_fibres=number_of_fibres, use_total=use_total

	version_module = '0.2'				; Version

	; If the version is requested, print it then quit
	if n_elements(version) ne 0 then begin
	    version = version_module
	    return
	endif

	;-------------------------------------------------------------------------------
	; Read the data ----------------------------------------------------------------

	; See the list of DRP extensions here:
	;	https://trac.sdss.org/wiki/MANGA/TRM/TRM_v1.0/datamodel

	; Relevant to this code are:
	;	1st extension is the flux
	data=double( READFITS(datacube_name,header,exten_no=1) )
	;	2nd extension is the inverse variance
	inv=double(readfits(datacube_name,exten_no=2))  
	error=sqrt(1./temporary(inv))			; Convert to error
							; TODO: Convert other code to use inv?
	;	4th extension is the wavelength solution
	wavelength=double(readfits(datacube_name,exten_no=4))
	;-------------------------------------------------------------------------------


	; Select wavelengths within the provided range
	lam_sel=indgen(n_elements(wavelength))
	if n_elements(lrange) ne 0 then begin
	    lam_sel = where(wavelength ge lrange[0] and wavelength le lrange[1])
	endif

	; Select RSS vs datacube based on the dimensionality of the data file
	;	s[0] = NAXIS ; s[j] = NAXISj
	sz=size(data)
	if sz[0] eq 3 then RSS = 0   ; datacube input
	if sz[0] eq 2 then RSS = 1   ; RSS input

	; TODO: Create separate modules for reading RSS vs. DATACUBES

	; if RSS, rearrange size (DRP output sets spectra along rows (not columns)
	if RSS eq 1 then begin
	    RESOLUTION_ELEMENT = 0.5			; set a fixed resolution element

	    data_ = transpose(temporary(data))		; flip the axes
	    error_ = transpose(temporary(error))

	    ; 7th and 8th extension are the x and y positions
	    x2d_=readfits(datacube_name,exten_no=7)
	    y2d_=readfits(datacube_name,exten_no=8)

	    x2d=median(temporary(x2d_[lam_sel,*]),/even,dimension=1)
	    y2d=median(temporary(y2d_[lam_sel,*]),/even,dimension=1)

	    ;--
	    ;   KEEPING FIBRES, DITHERS, AND EXPOSURES SEPARATED
;	    sz_=size(data_)
	    data  = temporary(data_)
	    error = temporary(error_)
	endif
	
	sz=size(data)

	; remove error spectra if bad for more than 20% of wavelength extension
	if RSS eq 0 then begin	; datacube file
	    for j = 0, sz[2]-1 do begin
		for i = 0 , sz[1] -1 do begin
		    indici = where(error[i,j,*] le 0 or finite(error[i,j,*]) ne 1 $
				   or finite(data[i,j,*]) ne 1)

		    if n_elements(indici) ge 0.2 * n_elements(data[i,j,*]) $
		       or min(data[i,j,*]) eq max(data[i,j,*]) then begin

;			data[i,j,*] = data[i,j,*]*0./0.
;			error[i,j,*] = error[i,j,*]*0.+1.
			data[i,j,*] = !VALUES.F_NAN
			error[i,j,*] = 1.0
		    endif
		endfor
	    endfor
	endif else begin	; RSS file
	    for i = 0,sz[1]-1 do begin
		indici = where(error[i,*] le 0 or finite(error[i,*]) ne 1 $
			       or finite(data[i,*]) ne 1)
		if n_elements(indici) ge 0.2 * n_elements(data[i,*]) $
		   or min(data[i,*]) eq max(data[i,*]) then begin

;		    data[i,*] = data[i,*]*0./0.
;		    error[i,*] = error[i,*]*0.+1.
		    data[i,*] = !VALUES.F_NAN
		    error[i,*] = 1.0
		endif
	    endfor
	endelse

	; Calculate S/N from data before resampling
	; CALLS MDAP_CALCULATE_SPECTRUM_SN for each spectrum
	if RSS eq 0 then begin      ;datacube
	    signal = fltarr(sz[1],sz[2])
	    noise = fltarr(sz[1],sz[2])
	    for j = 0, sz[2]-1 do begin
	        for i = 0, sz[1]-1 do begin
		    if ~keyword_set(use_total) then begin
		        MDAP_CALCULATE_SPECTRUM_SN, data[i,j,lam_Sel], error[i,j,lam_Sel], $
						    wavelength[lam_Sel], junk, signal=sss, $
						    noise=nnn
			signal[i,j]=sss[0]
			noise[i,j] = nnn[0]
		    endif else begin
			signal[i,j]=total(data[i,j,lam_Sel])
			noise[i,j] = sqrt(total(error[i,j,lam_Sel]^2))
		    endelse
		endfor
	    endfor
	endif else begin        ; rss
	    if ~keyword_set(use_total) then begin
		signal = median(data[*,lam_Sel],dimension=2,/even)
		noise  = median(error[*,lam_Sel],dimension=2,/even)
	    endif else begin
		nrss = n_elements(data[*,0])
		signal=fltarr(nrss)
		noise=fltarr(rss)
		for j = 0, nrss-1 do begin
		    signal[j] = total(data[j,lam_Sel])
		    noise[j] = sqrt(total(error[j,lam_Sel]^2))
		endfor
	    endelse 
	    
	    indici = where(finite(noise) eq 0 or finite(signal) eq 0  or noise le 0,compl=pos)
	    if indici[0] ne -1 then noise(indici) = max(noise[pos])
;	    if indici[0] ne -1 then signal(indici) = 0./0.
	    if indici[0] ne -1 then signal(indici) = !VALUES.F_NAN
	endelse


	; resample to a linear step in angstroms per pixel over a new wavelength range, if required
	if ~keyword_set(keep_original_step) then begin
	    min_step_lin = min( wavelength[1:n_elements(wavelength)-2] $
				- wavelength[0:n_elements(wavelength)-1])
	    nw=round((max(wavelength)-min(wavelength))/min_step_lin)
	    wavelength_linear_constant_step = dindgen(nw)*min_step_lin+min(wavelength)
	    
	    if RSS eq 0 then begin  ;datacube
	    
		data_interpolated = dblarr(sz[1],sz[2],n_elements(wavelength_linear_constant_step))
		error_interpolated = dblarr(sz[1],sz[2],n_elements(wavelength_linear_constant_step))
		
		for i = 0, sz[1]-1 do begin
		    for j = 0, sz[2]-1 do begin
		        data_interpolated[i,j,*] = interpol(data[i,j,*], wavelength, $
							    wavelength_linear_constant_step)
			error_interpolated[i,j,*] = interpol(error[i,j,*], wavelength, $
							     wavelength_linear_constant_step)
		    endfor
		endfor
		
		wavelength = temporary(wavelength_linear_constant_step)
		data = temporary(data_interpolated)
		error = temporary(error_interpolated)
	    endif else begin   ;RSS
		data_interpolated = dblarr(n_elements(wavelength_linear_constant_step),sz[2])
		error_interpolated = dblarr(n_elements(wavelength_linear_constant_step),sz[2])
		for i = 0, sz[1]-1 do begin
		    data_interpolated[i,*] = interpol(data[i,*], wavelength, $
						      wavelength_linear_constant_step)
		    error_interpolated[i,*] = interpol(error[i,*], wavelength, $
						       wavelength_linear_constant_step)
		endfor
		
		wavelength = temporary(wavelength_linear_constant_step)
		data = temporary(data_interpolated)
		error = temporary(error_interpolated)
	    endelse
	endif

	indx=where(error le 0 or finite(error) ne 1,compl=pos)
	
	if indx[0] ne -1 then error[indx] = median(error[pos])

	sz=size(data)

	; spatial coordinates
	if RSS eq 0 then begin	; datacubes
	    cdelt1=sxpar(header,'CDELT1')
	    cdelt2=sxpar(header,'CDELT2')
	    x=(findgen(sz[1])-(sz[1]-1)/2.)*cdelt1
	    y=(findgen(sz[2])-(sz[2]-1)/2.)*cdelt2
	    x2d=x#(y*0+1.)
	    y2d=(x*0.+1.)#y
	    x2d_reconstructed = x2d
	    y2d_reconstructed = y2d
	    signal2d_reconstructed = signal
	    noise2d_reconstructed = noise

	    header2d=["SIMPLE  =                    T / Written by IDL:",$
	              "BITPIX  =                  -32 / Number of bits per data pixel",$
		      "NAXIS   =                    3 / Number of data axes",$
		      "EXTEND  =                    T / FITS data may contain extensions",$
		      "CRPIX1  =              "+MDAP_STC((sz[1]+1)/2.)+" /",$
		      "CTYPE1  =                'arcsec' /",$
		      "CRVAL1  =                  0.0    /",$
		      "CDELT1  =          "+string(cdelt1)+" /",$
		      "CRPIX2  =              "+MDAP_STC((sz[2]+1)/2.)+" /",$
		      "CTYPE2  =                'arcsec' /",$
		      "CRVAL2  =                  0.0    /",$
		      "CDELT2  =          "+string(cdelt2)+" /",$
		      "OBJECT  =   "+datacube_name+"  /"]

	endif else begin	; for RSS, interpolate some info to 2D array
	    cdelt1 = RESOLUTION_ELEMENT
	    cdelt2 = RESOLUTION_ELEMENT
	    max_d = round(max([abs(x2d),abs(y2d)]))+2.
	    nbinsx = max_d*2/cdelt1
	    nbinsx = fix(nbinsx/2)*2+1
	    nbinsy = max_d*2/cdelt2
	    nbinsy = fix(nbinsy/2)*2+1

	    x=cdelt1*(findgen(nbinsx)-(nbinsx-1)/2)
	    y=cdelt2*(findgen(nbinsy)-(nbinsy-1)/2)

	    x2d_reconstructed=x#(y*0.+1.)
	    y2d_reconstructed=(x*0.+1.)#y
   
	    ; mdap_interpolate_2dmaps,signal,x2d,y2d,x2d_reconstructed,y2d_reconstructed,x2d_reconstructed,y2d_reconstructed,signal2d_reconstructed
	    ; signal2d_reconstructed = reform(signal2d_reconstructed,nbinsx,nbinsy)
	    ; mdap_interpolate_2dmaps,noise,x2d,y2d,x2d_reconstructed,y2d_reconstructed,x2d_reconstructed,y2d_reconstructed,noise2d_reconstructed
	    ; noise2d_reconstructed = reform(noise2d_reconstructed,nbinsx,nbinsy)
	    
;	    signal2d_reconstructed = x2d_reconstructed/0.
;	    noise2d_reconstructed =  x2d_reconstructed/0.
	    signal2d_reconstructed[*,*] = !VALUES.F_NAN
	    noise2d_reconstructed[*,*] = !VALUES.F_NAN
	    for j = 0, nbinsy-1 do begin
		for i = 0, nbinsx-1 do begin
		    dist = sqrt((x2d_reconstructed[i,j]-x2d)^2+(y2d_reconstructed[i,j]-y2d)^2)
		    min_dist = min(dist)
		    indici = where(dist le min_dist)
		    if min_dist le 1.5 then signal2d_reconstructed[i,j] = mean(signal[indici])
		    if min_dist le 1.5 then noise2d_reconstructed[i,j] = mean(noise[indici])
		endfor
	    endfor

	    header2d=["SIMPLE  =                    T / Written by IDL:",$
		      "BITPIX  =                  -32 / Number of bits per data pixel",$
		      "NAXIS   =                    3 / Number of data axes",$
		      "EXTEND  =                    T / FITS data may contain extensions",$
		      "CRPIX1  =              "+MDAP_STC((nbinsx+1)/2.)+" /",$
		      "CTYPE1  =                'arcsec' /",$
		      "CRVAL1  =                  0.0    /",$
		      "CDELT1  =          "+string(cdelt1)+" /",$
		      "CRPIX2  =              "+MDAP_STC((nbinsy+1)/2.)+" /",$
		      "CTYPE2  =                'arcsec' /",$
		      "CRVAL2  =                  0.0    /",$
		      "CDELT2  =          "+string(cdelt2)+" /",$
		      "OBJECT  =   "+datacube_name+"  /"]
	endelse

	;-- PATCH added from version 0.7 onwards, to center the coordinates
	;   onto the galaxy center. Header informations are updated accordingly.
	; get image center from signal2d_reconstructed
	sz_=size(signal2d_reconstructed)
	GCNTRD, signal2d_reconstructed, sz_[1]/2., sz_[2]/2., xcen, ycen, 3.
	xcen_pix=xcen
	ycen_pix=ycen
	x0=bilinear(x2d_reconstructed, xcen_pix, ycen_pix)
	xcen_pix=xcen
	ycen_pix=ycen
	y0=bilinear(y2d_reconstructed,xcen_pix,ycen_pix)
	if abs(x0) ge 4 then x0 = 0.
	if abs(y0) ge 4 then y0 = 0.

	x2d=x2d-x0
	y2d=y2d-y0
	x2d_reconstructed=x2d_reconstructed-x0
	y2d_reconstructed=y2d_reconstructed-y0
	header2d[6] = "CRVAL1  =                 "+MDAP_STC(-x0)+"    /"
	header2d[10] = "CRVAL2  =                "+MDAP_STC(-y0)+"    /"

end



