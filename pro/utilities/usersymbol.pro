;Brian Jackel   University of Western Ontario
;Bug reports cheerfully accepted
;jackel@canlon.physics.uwo.ca
;+
; NAME:         UserSymbol
;
; PURPOSE:      Make neat little user defined symbols
;
; CATEGORY:     Plotting/Graphics
;
; CALLING SEQUENCE:     UserSymbol,symbol_name
;
; INPUTS:
;       Symbol_name     a string containing the name of the desired symbol.
;                       Some possible options are Square, Triangle, Circle,
;                       Hexagon, BigX, Clover, Spiral, Star...
;
; KEYWORD PARAMETERS:
;       SIZE    Symbol size      (default=1)
;       LIST    if set, puts the list of available symbol names
;                 in the input parameter Symbol_Name
;       HELP    if set, returns this documentation header
;
; and also the keywords which apply to USERSYM
;       THICK   Line thickness   (default=1)
;       FILL    Fill symbol?     (default=0=no)
;       COLOR   Symbol color

;
; SIDE EFFECTS:         Calls  USERSYM to load the new symbol
;
; MODIFICATION HISTORY:      Brian Jackel   August 10 1992
;                            University of Western Ontario
;
; Bjj June 2 1994      Fixed up the handling of no clear match.
;-
pro USERSYMBOL,symbol_name,SIZE_OF_SYMBOL=size_of_symbol, $
                           ORIENTATION=orientation,       $
                           LIST=list,HELP=help,_EXTRA=_extra

  ON_ERROR,2

  IF KEYWORD_SET(HELP) THEN BEGIN
     DOC_LIBRARY,'USERSYMBOL'
     RETURN
  ENDIF

  symbol_list= ['DIAMOND','PENTAGON','CLOVER','PACMAN','SPIRAL','BIGX']
  symbol_list= [symbol_list,'CIRCLE','SQUARE','TRIANGLE','STAR','HEXAGON']

  IF KEYWORD_SET(LIST) THEN BEGIN
      symbol_name= symbol_list
      return                   ;return a list of the available symbols
  ENDIF

  IF not KEYWORD_SET(SIZE_OF_SYMBOL) THEN symsize=!p.symsize ELSE symsize= (size_of_symbol > 0.01) < 100.0
  IF (symsize EQ 0) THEN symsize= 1.0   ;because !p.symsize is sometimes zero

  symbol= STRUPCASE( STRCOMPRESS(symbol_name,/REMOVE_ALL) )

  CASE symbol OF
  'DIAMOND': BEGIN
               x= [0.0,0.8,0.0,-0.8,0.0]
               y= [1.2,0.0,-1.2,0.0,1.2]
             END
  'PENTAGON':BEGIN
               theta= findgen(6)/5 * 360.0 * !dtor
               x= sin(theta)
               y= cos(theta)
             END
  'CLOVER':  BEGIN
               theta= findgen(41)/40.0 * 360.0 * !dtor
               r= ABS(1.0 *symsize* sin(2.0*theta))
               x= r * sin(theta)
               y= r * cos(theta)
             END
  'PACMAN':  BEGIN
               theta= (- findgen(41)/50.0*360.0 + 35.0 )*!dtor
               x= [0.0, sin(theta), 0.0]
               y= [0.0, cos(theta) ,0.0]
             END
  'SPIRAL':  BEGIN
               theta= findgen(41)/40.0 * 720.0 * !dtor
               r= theta / MAX(theta)
               x= r * sin(theta)
               y= r * cos(theta)
             END
  'BIGX':    BEGIN
               x= 0.34 * [0,2,3,3,1, 3, 3, 2, 0,-2,-3,-3,-1,-3,-3,-2,0]
               y= 0.34 * [1,3,3,2,0,-2,-3,-3,-1,-3,-3,-2, 0, 2, 3, 3,1]
             END
  'CIRCLE':  BEGIN
               n= 17.0
               theta= findgen(n)/(n-1.0) * 360.0 * !dtor
               x= sin(theta)
               y= cos(theta)
             END
  'SQUARE':  BEGIN
               theta= (findgen(5)/4.0 * 360.0 + 45.0 )*!dtor
               x= sin(theta)
               y= cos(theta)
             END
  'TRIANGLE':BEGIN
               theta= [0,120,240,360]*!dtor
               x= sin(theta)
               y= cos(theta)
             END
  'STAR':    BEGIN
               theta= [0,36, 72,108, 144,180, 216,252, 288,324,0]*!dtor
               r= [1.0,0.4, 1.0,0.4, 1.0,0.4, 1.0,0.4, 1.0,0.4,1.0]
               x= r *sin(theta)
               y= r *cos(theta)
             END
  'HEXAGON': BEGIN
               theta= [0,60,120,180,240,300,360]*!dtor
               x= sin(theta)
               y= cos(theta)
             END
  'SPIRAL2': BEGIN
              n=49
              theta= 2.0*!pi*FINDGEN(n)/((n-1)/2.0)
              r= FINDGEN(n)/(n-1)
              x= r*SIN(theta)
              y= r*COS(theta)
            END
  ELSE: BEGIN
         MESSAGE,'Unrecognized symbol name, searching for match',/INFORMATIONAL
         hits= STRPOS( symbol_list, symbol )
         w= WHERE(hits NE -1)
         IF (w(0) NE -1) THEN BEGIN ;at least one substring match, use
            hit_names= symbol_list(w(0))
            FOR i=1,n_elements(w)-1 DO hit_names= hit_names +   $
                                                  ' ' + symbol_list(w(i))
            MESSAGE,'...possible matches: '+hit_names,/INFORMATIONAL
            MESSAGE,'...will use the first (or only) one',/INFORMATIONAL
            symbol_name= symbol_list(w(0))  ;recursion to help us out
            USERSYMBOL,symbol_name,_EXTRA=_extra
         ENDIF ELSE BEGIN
            MESSAGE,'...no quick match.  Try USERSYMBOL,list,/LIST',/INFORMATIONAL
         ENDELSE
         return   ;either with a guessed symbol, or a list of them
        END
  ENDCASE


;Introduce scaling to the symbol size, if requested
;
IF (symsize NE 1.0) THEN BEGIN
   x= x * symsize
   y= y * symsize
ENDIF

;Rotate the symbol, if requested
;
 IF KEYWORD_SET(ORIENTATION) THEN BEGIN
    r= SQRT(x^2 + y^2)
    theta= ATAN(y,x)
    theta= theta + orientation*!dtor
    x= r * COS(theta)
    y= r * SIN(theta)
 ENDIF


;Use the library routine USERSYM to set up the symbol
;
  USERSYM,x,y,_EXTRA=_extra

  RETURN
END
