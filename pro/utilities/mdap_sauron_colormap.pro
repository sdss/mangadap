;----------------------------------------------------------------------
pro mdap_sauron_colormap
;
; The official SAURON colormap: type "sauron_colormap" to load the colormap.
;
; Michele Cappellari & Eric Emsellem, Leiden, 10 July 2001
;
; Start with these 7 equally spaced coordinates, then add 4 additional points
; x = findgen(7)*255/6. + 1
; 1.0  43.5  86.0  128.5  171.0  213.5  256.0
;
x = [1.0, 43.5, 86.0, 86.0+20, 128.5-10, 128.5, 128.5+10, 171.0-20, 171.0, 213.5, 256.0]

red =   [0.0, 0.0, 0.4,  0.5, 0.3, 0.0, 0.7, 1.0, 1.0,  1.0, 0.9]
green = [0.0, 0.0, 0.85, 1.0, 1.0, 0.9, 1.0, 1.0, 0.85, 0.0, 0.9]
blue =  [0.0, 1.0, 1.0,  1.0, 0.7, 0.0, 0.0, 0.0, 0.0,  0.0, 0.9]

xnew = findgen(256)+1

redV = INTERPOL(red, x, xnew)
greenV = INTERPOL(green, x, xnew)
blueV = INTERPOL(blue, x, xnew)

tvlct, redV*255, greenV*255, blueV*255 ; load the SAURON colormap

end
;----------------------------------------------------------------------

