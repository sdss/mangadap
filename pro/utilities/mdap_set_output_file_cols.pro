;+
; NAME:
;       MDAP_SET_DRPS_COLS
;       MDAP_SET_BINS_COLS
;       MDAP_SET_ELPAR_COLS
;       MDAP_SET_STFIT_COLS
;       MDAP_SET_SGFIT_COLS
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
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Functions that declare the name of the column names of the binary tables
FUNCTION MDAP_SET_DRPS_COLS
        return, [ 'XPOS', 'YPOS', 'SIGNAL', 'NOISE', 'BINID', 'BINW' ]
END
FUNCTION MDAP_SET_BINS_COLS
        return, [ 'XBIN', 'YBIN', 'BINA', 'BINSN', 'NBIN', 'BINF' ]
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
FUNCTION MDAP_SET_SIPAR_COLS
        return, [ 'SINAME', 'PASSBAND', 'BLUEBAND', 'REDBAND', 'UNIT' ]
END
FUNCTION MDAP_SET_SINDX_COLS
        return, [ 'SIOMIT', 'INDX', 'INDXERR', 'INDX_OTPL', 'INDX_BOTPL' ]
END

