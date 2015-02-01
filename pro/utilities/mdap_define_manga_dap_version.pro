;+
; NAME:
;       MDAP_DEFINE_MANGA_DAP_VERSION
;
; PURPOSE:
;       Define the structure that holds all the version infomation for
;       the DAP.
;
; CALLING SEQUENCE:
;       result = MDAP_DEFINE_MANGA_DAP_VERSION()
;
; OUTPUT:
;       result is the default version structure, where all the versions
;       are set to '0'.
;
; REVISION HISTORY:
;       30 Jan 2015: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_DEFINE_MANGA_DAP_VERSION

        version = { MaNGADAPVersion, main:'0', execute_plan:'0', read_drp_fits:'0', $
                                     fiducial_bin_xy:'0', tpl_lib_setup:'0', calculate_sn:'0', $
                                     spatial_binning:'0', spectral_fitting:'0', $
                                     emission_line:'0', spectral_index:'0' }

        return, version
END

