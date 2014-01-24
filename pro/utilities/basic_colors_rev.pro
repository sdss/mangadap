;#> basic_colors.dc1
; Identifier   basic_colors
;
; Purpose      define basic colors
;
; Synopsis     basic_colors [, black , white , red    , green    , blue $
;                            , yellow, cyan  , magenta, orange   , mint $
;                            , purple, pink  , olive  , lightblue, gray $
;                            , /test ]   
;
; Arguments    Name   I/O Type:     Description:
;	       ----------------------------------------------------------
;	       black   O  int	    color index black
;	       ...     O  int	    ...
;              /test   I            make an example plot of all colours
;
; Description  defines B&W, basic RGB colors, lightb(lue) and gray
;
; Returns      literal keywords for color indices

;basic_colors_rev, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   

; Example      basic_colors_rev,   black , white , red    , green    , blue, $
;                              yellow, cyan  , magenta, orange   , mint, $
;                              purple, pink  , olive  , lightblue, gray   
;
;	       plot ,[0,1],[0,4],color=yellow
;	       oplot,[0,1],[0,3],color=blue
;	       oplot,[0,1],[0,2],color=green
;	       oplot,[0,1],[0,1],color=red
;
; Category     UTIL
;
; Filename     basic_colors.pro
;
; Author       D. Kunze (DIK)
;
; Version      1.5
;
; History      1.0 03-04-98 DIK   insert IA module
;              1.1 10-10-00 FL/DK SPR_S0579 update colors
;              1.2 10-10-00 FL    bottom parameter to loadct
;              1.3 23-01-01 FL    only load color table once
;              1.4 18-04-02 FL    true colors and added test option
;              1.5 12-08-02 FL    decomposed=0 only for X device
;              
; -----------------------------------------------------------------------------
; Copyright (C) 1998, Max-Planck-Institut fuer extraterrestrische Physik (MPE)
; Garching, Germany. All rights reserved. Unauthorized reproduction prohibited.
; -----------------------------------------------------------------------------
;#<
PRO basic_colors_rev,   black , white , red    , green , blue, $
                    yellow, cyan  , magenta, orange, mint, $
                    purple, pink  , olive  , lightb, gray, $
                    test = test

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr   ;for looadct

;	 0 1 2 3 4 5 6 7  8   9   10  11  12  13  14

; set this to get pseudo colors on a true color device
; on a pseudo color device this keyword has no effect
;
  if !d.name EQ 'X' then device,decomposed=0
;
; load a color table if this has not yet been done
; keep first fifteen colors free
;
 
   if n_elements(r_curr) le 0 then loadct,39,bottom=15

 red  =[1,0,1,0,0,1,0,1,1.0,0.3,0.7,1.0,0.7,0.3,0.5]
 green=[1,0,0,1,0,1,1,0,0.7,1.0,0.3,0.3,1.0,0.7,0.5]
 blue =[1,0,0,0,1,0,1,1,0.3,0.7,1.0,0.7,0.3,1.0,0.5]
;   red  =[0,1,1,0,0,1,0,1,1.0,0.3,0.7,1.0,0.7,0.3,0.5]
;   green=[0,1,0,1,0,1,1,0,0.7,1.0,0.3,0.3,1.0,0.7,0.5]
;   blue =[0,1,0,0,1,0,1,1,0.3,0.7,1.0,0.7,0.3,1.0,0.5]

  tvlct,255*red,255*green,255*blue

; NB FL, colours as they 'should' appear on a TrueColor device
; black, white, red, green, blue, yellow, cyan, magenta, $
; orange, mint, purple, pink, olive, lightb, grey
;
  colours = [ 'white' , 'black' , 'red'    , 'green' , 'blue', $
              'yellow', 'cyan'  , 'magenta', 'orange', 'mint', $
              'purple', 'pink'  , 'olive'  , 'lightb', 'gray' ]

  black  = 1
  white  = 0
  red	 = 2
  green  = 3
  blue	 = 4
  yellow = 5
  cyan	 = 6
  magenta= 7
  orange = 8
  mint	 = 9
  purple = 10
  pink   = 11
  olive  = 12
  lightb = 13
  gray	 = 14

  if keyword_set(test) then begin
    if !d.name eq 'X' then $
         plot,indgen(10)+1,psym=0,xrange=[0,12],yrange=[0,25], $
              color=1,thick=2.5 $
    else plot,indgen(10)+1,psym=0,xrange=[0,12],yrange=[0,25], $
              color=0,thick=2.5
    for col=0,14 do begin
      oplot,indgen(10)+col,psym=0,color=col,thick=2.5
      xyouts,9.5,col+9,colours[col],/data,color=col,charsize=1.25
    endfor
  endif

  RETURN

END
