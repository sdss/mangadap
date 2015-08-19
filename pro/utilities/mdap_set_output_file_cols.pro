;+
; NAME:
;       MDAP_SET_DRPS_COLS
;       MDAP_SET_BINS_COLS
;       MDAP_SET_ELPAR_COLS
;       MDAP_SET_STFIT_COLS
;       MDAP_SET_SGFIT_COLS
;       MDAP_SET_ELOPAR_COLS
;       MDAP_SET_ELOFIT_COLS
;       MDAP_SET_SIPAR_COLS
;       MDAP_SET_SINDX_COLS
;
; PURPOSE:
;       A set of functions that define the different sets of binary table
;       columns names.  Include these functions in a procedure file using:
;
;       @mdap_set_output_file_cols
;
; REVISION HISTORY:
;       07 Nov 2014: Original implementation by K. Westfall (KBW)
;       06 Dec 2014: (KBW) Changes to DRPS (includes BINVR)
;       12 Dec 2014: (KBW) Add columns for ELOPAR and ELOFIT tables
;       09 Jan 2015: (KBW) Include instrumental dispersion for GANDALF fit
;       09 Feb 2015: (KBW) Include fraction of good pixels and min(flux)
;                          == max(flux) flag in DRPS extension.
;       13 Jul 2015: (KBW) Changes to emission-line-only output columns
;       12 Aug 2015: (KBW) Added new columns for spectral index output.
;       18 Aug 2015: (KBW) Added functions for the non-parametric
;                          emission-line measurement columns.
;-
;------------------------------------------------------------------------------

; !! IF THESE CHANGE, ALSO NEED TO CHECK IF MDAP_ANALYSIS_BLOCK_LOGIC
; !! NEEDS TO CHANGE

;-------------------------------------------------------------------------------
; Functions that declare the name of the column names of the binary tables
FUNCTION MDAP_SET_DRPS_COLS
        return, [ 'XPOS', 'YPOS', 'FGOODPIX', 'MINEQMAX', 'SIGNAL', 'NOISE', 'BINVR', 'BINID', $
                  'BINW' ]
END
FUNCTION MDAP_SET_BINS_COLS
        return, [ 'BINXRL', 'BINYRU', 'BINR', 'BINA', 'BINSN', 'NBIN', 'BINF' ]
END
FUNCTION MDAP_SET_ELPAR_COLS
        return, [ 'ELNAME', 'RESTWAVE', 'TIEDKIN', 'TIEDTYPE', 'DOUBLET' ]
END
FUNCTION MDAP_SET_STFIT_COLS
        return, [ 'TPLW', 'ADDPOLY', 'MULTPOLY', 'KIN', 'KINERR', 'RCHI2' ]
END
FUNCTION MDAP_SET_SGFIT_COLS
        return, [ 'TPLW', 'MULTPOLY', 'KIN', 'KINERR', 'RCHI2', 'RED', 'REDERR', 'ELOMIT', 'AMPL', $
                  'AMPLERR', 'IKIN', 'IKINERR', 'SINST', 'FLUX', 'FLUXERR', 'EW', 'EWERR' ]
END
FUNCTION MDAP_SET_ELBAND_COLS
        return, [ 'ELNAME', 'RESTWAVE', 'BANDPASS', 'BLUESIDE', 'REDSIDE' ]
END
FUNCTION MDAP_SET_ELMMNT_COLS
        return, [ 'OMIT', 'FLUX_RAW', 'FLUXERR_RAW', 'MOM1_RAW', 'MOM1ERR_RAW', 'MOM2_RAW', $
                  'MOM2ERR_RAW', 'BLUF_RAW', 'BLUFERR_RAW', 'REDF_RAW', 'REDFERR_RAW', 'BFCEN', $
                  'BCONT', 'BCONTERR', 'RFCEN', 'RCONT', 'RCONTERR', 'FLUX', 'FLUXERR', 'MOM1', $
                  'MOM1ERR', 'MOM2', 'MOM2ERR' ]
END
FUNCTION MDAP_SET_ELOPAR_COLS
        return, [ 'ELNAME', 'RESTWAVE' ]
END
FUNCTION MDAP_SET_ELOFIT_COLS
        return, [ 'KIN_EW', 'KINERR_EW', 'KINSTDE_EW', 'NKIN_EW', 'ELOMIT_EW', 'AMPL_EW', $
                  'AMPLERR_EW', 'IKIN_EW', 'IKINERR_EW', 'SINST_EW', 'FLUX_EW', 'FLUXERR_EW', $
                  'EW_EW', 'EWERR_EW', $
                  'KIN_FB', 'KINERR_FB', 'KINSTDE_FB', 'NKIN_FB', 'ELOMIT_FB', 'AMPL_FB', $
                  'AMPLERR_FB', 'IKIN_FB', 'IKINERR_FB', 'SINST_FB', 'FLUX_FB', 'FLUXERR_FB', $
                  'FB_FB', 'FBERR_FB']
END
FUNCTION MDAP_SET_SIPAR_COLS
        return, [ 'SINAME', 'PASSBAND', 'BLUEBAND', 'REDBAND', 'UNIT' ]
END
FUNCTION MDAP_SET_SINDX_COLS
        return, [ 'SIOMIT', 'BCONT', 'BCONTERR', 'RCONT', 'RCONTERR', 'INDX_RAW', 'INDXERR_RAW', $
                  'INDXPERR_RAW', 'BCONT_OTPL', 'RCONT_OTPL', 'INDX_OTPL', 'BCONT_BOTPL', $
                  'RCONT_BOTPL', 'INDX_BOTPL', 'INDX', 'INDXERR', 'INDXPERR' ]
        
END

