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
;       07 Nov 2014: (KBW) Original implementation
;       06 Dec 2014: (KBW) Changes to DRPS (includes BINVR)
;       12 Dec 2014: (KBW) Add columns for ELOPAR and ELOFIT tables
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Functions that declare the name of the column names of the binary tables
FUNCTION MDAP_SET_DRPS_COLS
        return, [ 'XPOS', 'YPOS', 'SIGNAL', 'NOISE', 'BINVR', 'BINID', 'BINW' ]
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
                  'AMPLERR', 'IKIN', 'IKINERR', 'FLUX', 'FLUXERR', 'EW', 'EWERR' ]
END
FUNCTION MDAP_SET_ELOPAR_COLS
        return, [ 'ELNAME', 'RESTWAVE' ]
END
FUNCTION MDAP_SET_ELOFIT_COLS
        return, [ 'KIN_EW', 'KINERR_EW', 'ELOMIT_EW', 'AMPL_EW', 'AMPLERR_EW', 'IKIN_EW', $
                  'IKINERR_EW', 'SINST_EW', 'FLUX_EW', 'FLUXERR_EW', 'EW_EW', 'EWERR_EW', $
                  'KIN_FB', 'KINERR_FB', 'ELOMIT_FB', 'AMPL_FB', 'AMPLERR_FB', 'IKIN_FB', $
                  'IKINERR_FB', 'SINST_FB', 'FLUX_FB', 'FLUXERR_FB', 'EW_FB', 'EWERR_FB' ]
END
FUNCTION MDAP_SET_SIPAR_COLS
        return, [ 'SINAME', 'PASSBAND', 'BLUEBAND', 'REDBAND', 'UNIT' ]
END
FUNCTION MDAP_SET_SINDX_COLS
        return, [ 'SIOMIT', 'INDX', 'INDXERR', 'INDX_OTPL', 'INDX_BOTPL' ]
END

