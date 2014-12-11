;+
; NAME:
;       MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE
;
; PURPOSE:
;       Define the hard-coded wavelengths used by
;       MDA_EMISSION_LINE_ONLY_FIT to fit the emission-line only
;       spectrum using either Enci's or Francesco's code.
;
; CALLING SEQUENCE:
;       result = MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()
;
; OUTPUT:
;       The result is an array of EmissionLine structures.  See
;       MDAP_ASSIGN_EMISSION_LINE_PARAMETERS.
;
; PROCEDURES CALLED:
;       MDAP_ASSIGN_EMISSION_LINE_PARAMETERS()
;
; REVISION HISTORY:
;       11 Dec 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------


FUNCTION MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE

        ; The fitted lines and their vacuum rest wavelengths:
        OII_3727_rest_wave = 3726.03
        AIRTOVAC, OII_3727_rest_wave

        Hb_4861_rest_wave = 4861.32
        AIRTOVAC, Hb_4861_rest_wave

        OIII_4959_rest_wave = 4958.83
        AIRTOVAC, OIII_4959_rest_wave

        OIII_5007_rest_wave = 5006.77
        AIRTOVAC, OIII_5007_rest_wave

        NII_6548_rest_wave = 6547.96
        AIRTOVAC, NII_6548_rest_wave

        Ha_6563_rest_wave = 6562.80
        AIRTOVAC, Ha_6563_rest_wave

        NII_6583_rest_wave = 6583.34
        AIRTOVAC, NII_6583_rest_wave

        SII_6717_rest_wave = 6716.31
        AIRTOVAC, SII_6717_rest_wave

        SII_6730_rest_wave = 6730.68
        AIRTOVAC, SII_6730_rest_wave

        ; Create the eml_par equivalent for these lines
        id = indgen(9)
        name = [ 'OII', 'Hb', 'OIII', 'OIII', 'NII', 'Ha', 'NII', 'SII', 'SII' ]
        lambda = [ OII_3727_rest_wave, Hb_4861_rest_wave, OIII_4959_rest_wave, NII_6548_rest_wave, $
                   Ha_6563_rest_wave, NII_6583_rest_wave, SII_6717_rest_wave, SII_6730_rest_wave ]
        action = [ 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f' ]
        kind = [ 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l' ]
        a = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]
        v = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        s = [ 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 ]
        fit = [ 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f' ]

        return, MDAP_ASSIGN_EMISSION_LINE_PARAMETERS(id, name, lambda, action, kind, a, v, s, fit)
END


