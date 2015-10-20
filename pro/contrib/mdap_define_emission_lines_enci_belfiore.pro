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
;       11 Dec 2014: Original implementation by K. Westfall (KBW)
;       13 Aug 2015: (KBW) Include second OII line, and the OI doublet.
;       20 Oct 2015: (KBW) Include a single line fit to the OII doublet.
;-
;------------------------------------------------------------------------------


FUNCTION MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE

        ; The fitted lines and their vacuum rest wavelengths:
        OII_3727_rest_wave = 3727.092d
        OII_3729_rest_wave = 3729.875d
        OII_3727_doublet_rest_wave = (OII_3727_rest_wave + OII_3729_rest_wave)/2.
        Hb_4861_rest_wave = 4862.691d
        OIII_4959_rest_wave = 4960.295d
        OIII_5007_rest_wave = 5008.240d
        OI_6300_rest_wave = 6302.046d
        OI_6363_rest_wave = 6365.535d
        NII_6548_rest_wave = 6549.86d
        Ha_6563_rest_wave = 6564.632d
        NII_6583_rest_wave = 6585.271d
        SII_6717_rest_wave = 6718.294d
        SII_6730_rest_wave = 6732.674d

        ; Create the eml_par equivalent for these lines
        neml = 13
        id = indgen(neml)
        name = [ 'OII', 'OII', 'OIId', 'Hb', 'OIII', 'OIII', 'OI', 'OI', 'NII', 'Ha', 'NII', $
                 'SII', 'SII' ]
        lambda = [ OII_3727_rest_wave, OII_3729_rest_wave, OII_3727_doublet_rest_wave, $
                   Hb_4861_rest_wave, OIII_4959_rest_wave, OIII_5007_rest_wave, OI_6300_rest_wave, $
                   OI_6363_rest_wave, NII_6548_rest_wave, Ha_6563_rest_wave, NII_6583_rest_wave, $
                   SII_6717_rest_wave, SII_6730_rest_wave ]
        action = [ 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f' ]
        kind = [ 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l', 'l' ]
        a = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]
        v = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        s = [ 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 ]
        fit = [ 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f' ]

        return, MDAP_ASSIGN_EMISSION_LINE_PARAMETERS(id, name, lambda, action, kind, a, v, s, fit)
END


