;+
; NAME:
;       MDAP_DEFINE_NONPAR_EMISSION_LINE_BANDS
;
; PURPOSE:
;       Define the hard-coded wavelengths used by
;       MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS to perform non-parametric
;       measurements of the emission-lines.
;
;       WAVELENGTHS ARE DEFINED IN VACUUM!
;
; CALLING SEQUENCE:
;       result = MDAP_DEFINE_NONPAR_EMISSION_LINE_BANDS()
;
; OUTPUT:
;       The result is an array of EmissionLineBand structures.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;       17 Aug 2015: Original implementation by K. Westfall (KBW)
;       18 Aug 2015: (KBW) Changed structure element names to match
;                          EmissionLine structure; eases writting the
;                          band data do the output file.
;-
;------------------------------------------------------------------------------
FUNCTION MDAP_DEFINE_NONPAR_EMISSION_LINE_BANDS

        ; Only 11 bands because OII 3727 and 3729 are combined into a
        ;   single band

        ; Declare and replicate the structure
        neml = 11
        eml_bands = replicate( { EmissionLineBand, name:'', lambda:0.0d, bandpass:dblarr(2), $
                                                   blueside:dblarr(2), redside:dblarr(2) }, neml )

        j = 0

        ; OII doublet is a single bandpass
        eml_bands[j].name = 'OII'
        eml_bands[j].lambda = 3727.092d
        eml_bands[j].bandpass = [ 3716.3, 3738.3 ]
        eml_bands[j].blueside = [ 3696.3, 3716.3 ]
        eml_bands[j].redside  = [ 3738.3, 3758.3 ]
        j = j+1

        eml_bands[j].name = 'Hb'
        eml_bands[j].lambda = 4862.691d
        eml_bands[j].bandpass = [ 4852.7, 4872.7 ]
        eml_bands[j].blueside = [ 4798.9, 4838.9 ]
        eml_bands[j].redside  = [ 4885.6, 4925.6 ]
        j = j+1

        eml_bands[j].name = 'OIII'
        eml_bands[j].lambda = 4960.295d
        eml_bands[j].bandpass = [ 4950.3, 4970.3 ]
        eml_bands[j].blueside = [ 4930.3, 4950.3 ]
        eml_bands[j].redside  = [ 4970.3, 4990.3 ]
        j = j+1

        eml_bands[j].name = 'OIII'
        eml_bands[j].lambda = 5008.240d
        eml_bands[j].bandpass = [ 4998.2, 5018.2 ]
        eml_bands[j].blueside = [ 4978.2, 4998.2 ]
        eml_bands[j].redside  = [ 5018.2, 5038.2 ]
        j = j+1

        eml_bands[j].name = 'OI'
        eml_bands[j].lambda = 6302.046d
        eml_bands[j].bandpass = [ 6292.0, 6312.0 ]
        eml_bands[j].blueside = [ 6272.0, 6292.0 ]
        eml_bands[j].redside  = [ 6312.0, 6332.0 ]
        j = j+1

        eml_bands[j].name = 'OI'
        eml_bands[j].lambda = 6365.535d
        eml_bands[j].bandpass = [ 6355.5, 6375.5 ]
        eml_bands[j].blueside = [ 6335.5, 6355.5 ]
        eml_bands[j].redside  = [ 6375.5, 6395.5 ]
        j = j+1

        eml_bands[j].name = 'NII'
        eml_bands[j].lambda = 6549.86d
        eml_bands[j].bandpass = [ 6542.9, 6556.9 ]
        eml_bands[j].blueside = [ 6483.0, 6513.0 ]
        eml_bands[j].redside  = [ 6623.0, 6653.0 ]
        j = j+1

        eml_bands[j].name = 'Ha'
        eml_bands[j].lambda = 6564.632d
        eml_bands[j].bandpass = [ 6557.6, 6571.6 ]
        eml_bands[j].blueside = [ 6483.0, 6513.0 ]
        eml_bands[j].redside  = [ 6623.0, 6653.0 ]
        j = j+1

        eml_bands[j].name = 'NII'
        eml_bands[j].lambda = 6585.271d
        eml_bands[j].bandpass = [ 6575.3, 6595.3 ]
        eml_bands[j].blueside = [ 6483.0, 6513.0 ]
        eml_bands[j].redside  = [ 6623.0, 6653.0 ]
        j = j+1

        eml_bands[j].name = 'SII'
        eml_bands[j].lambda = 6718.294d
        eml_bands[j].bandpass = [ 6711.3, 6725.3 ]
        eml_bands[j].blueside = [ 6673.0, 6703.0 ]
        eml_bands[j].redside  = [ 6748.0, 6778.0 ]
        j = j+1

        eml_bands[j].name = 'SII'
        eml_bands[j].lambda = 6732.674d
        eml_bands[j].bandpass = [ 6725.7, 6739.7 ]
        eml_bands[j].blueside = [ 6673.0, 6703.0 ]
        eml_bands[j].redside  = [ 6748.0, 6778.0 ]
        j = j+1

        return, eml_bands
END


