;+
; NAME:
;       MDAP_ALL_SPECTRA_BINNING
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_ALL_SPECTRA_BINNING, xcoo, ycoo, signal, noise, binned_indx, binned_xcoo, $
;                                 binned_ycoo, binned_ston, nbinned, $
;                                 noise_calib=noise_calib, optimal_weighting=optimal_weighting
;
; INPUTS:
;       xcoo dblarr[N]
;               On-sky X coordinate for each of the N spectra.
;
;       ycoo dblarr[N]
;               On-sky Y coordinate for each of the N spectra.
;
;       signal dblarr[N]
;               Estimate of the spectrum signal.
;        
;       noise dblarr[N]
;               Estimate of the spectrum noise.
;
; OPTIONAL INPUTS:
;       noise_calib int
;               Flag to calibrate the noise vector to nominally account
;               for covariance.
;
;       optimal_weighting
;               Flag used to set optimal weighting.  If it exists, the
;               spectra are expected to be combined using S/(N)^2
;               weighting.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       binned_indx intarr[N]
;               Indicates in which bin each of the N spectra were
;               placed; all spectra are placed in the same bin, i=0.
;
;       binned_xcoo dblarr[1]
;               Luminosity-weighted on-sky X coordinate.
;
;       binned_ycoo dblarr[1]
;               Luminosity-weighted on-sky Y coordinate.
;
;       binned_ston dblarr[1]
;               S/N of the combined spectrum.
;
;       nbinned intarr[1]
;               Number of binned spectra.
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
;       08 Dec 2014: Original implementation by K. Westfall (KBW)
;       16 Mar 2015: (KBW) Remove sn_calibration vector in favor of
;                          noise_calib.
;-
;------------------------------------------------------------------------------

PRO MDAP_ALL_SPECTRA_BINNING, xcoo, ycoo, signal, noise, binned_indx, binned_xcoo, binned_ycoo, $
                              binned_ston, nbinned, noise_calib=noise_calib, $
                              optimal_weighting=optimal_weighting

        ; Set all spectra to be in the same bin
        ns = n_elements(xcoo)                       ; Number of spectra
        binned_indx = lonarr(ns)                    ; All spectra in bin 0

        wsum = total(signal)                        ; Calculate the luminosity-weighted coordinates
        binned_xcoo = [ total(signal*xcoo)/wsum ]
        binned_ycoo = [ total(signal*ycoo)/wsum ]
       
        ; Calculate the S/N
        binned_ston = [ MDAP_CALCULATE_BIN_SN(signal, noise, noise_calib=noise_calib, $
                                              optimal_weighting=optimal_weighting) ]

        nbinned = [ns]              ; Save the number of binned spectra
END

