;+
; NAME:
;       MDAP_INSTR_DISPERSION_AT_EMISSION_LINE
;
; PURPOSE:
;       Interpolate the instrumental dispersion (sigma) at the fitted
;       centroid of the emission lines.
;
; CALLING SEQUENCE:
;       MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, eml_par, omitted, kinematics, sinst
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength at each of the C wavelength channels.
;
;       sres dblarr[C]
;               Spectral resolution (R) at each wavelength channel.
;
;       eml_par EmissionLine[E]
;               Array of emission-line structures used in fitting the
;               emission lines.
;
;       omitted intarr[B][E]
;               Omission flags for the E emission lines when fitted to
;               each of the B binned spectra.
;
;       kinematics dblarr[B][E][K]
;               The K fitted kinematic moments for each of the E
;               emission lines in the B binned spectra.  The ONLY
;               element used is the velocity which MUST be at index K=0.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       sinst dblarr[B][E]
;               Instrumental dispersion (sigma in km/s) at the centroids
;               of the E emission lines fitted in the B binned spectra.
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
;       15 Dec 2014: (KBW) Original implementation
;       09 Jan 2015: (KBW) Use new MDAP_INSTRUMENTAL_DISPERSION function
;       10 Jul 2015: (KBW) Change default instrumental dispersion to
;                          -9999.0d
;-
;------------------------------------------------------------------------------

PRO MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, $
                wave, sres, eml_par, omitted, kinematics, sinst

        sz = size(kinematics)                       ; Dimensions of the kinematics array
        nbin = sz[1]                                ; Number of binned spectra
        neml = n_elements(eml_par)                  ; Number of emission lines
        if sz[2] ne neml then $
            message, 'Incorrect number of kinematic measurements entered!'

;       sinst = dblarr(nbin, neml)                  ; Initialize the instrumental dispersion
        ; Initialize the instrumental dispersion
        sinst = make_array(nbin, neml, /double, value=-9999.0d)
        for i=0,nbin-1 do begin
            indx = where(omitted[i,*] lt 1, count)
            if count eq 0 then $
                continue
            sinst[i,indx] = MDAP_INSTRUMENTAL_DISPERSION(wave, sres, eml_par[indx].lambda, $
                                                         kinematics[i,indx,0])
        endfor
END


