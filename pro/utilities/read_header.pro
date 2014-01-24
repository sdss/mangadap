pro read_header,filename,campo,valore,first_occur=first_occur,silent=silent,from_header=from_header
;filename: [string] name of the fits file to read
;campo:    [string] field to look for in the header
;valore:   [string] If campo is found, valore is the value of that field.
;                   If campo is not found, valore = -1
;first_occur  if set, the "valore" value of the first "campo" field is
;           read (in the case there are multiple entries for
;           "campo"). The default is to read the last "campo" field
;           value.
 
;from_header  the input header; if specified it overrides the header from
;             filename. default read the header from filename

;ima=readfits(filename,header,/silent)
if n_elements(from_header) eq 0 then begin
   header=headfits(filename)
endif else begin
   header=from_header
endelse

record_name_length=STRLEN(campo)
valore='-1'
valore_all='-1'
FOR i = 0, n_elements(header)-1 DO BEGIN
    key=STRMID(header[i],0,record_name_length)
    IF (key EQ campo) THEN BEGIN
        valore_all=header[i]
    if keyword_set(first_occur) then goto, esci_loop
     ENDIF
   
ENDFOR
esci_loop:
if valore_all eq '-1' then begin
 if not keyword_set(silent) then  print, 'field ',campo,' not found.  returning -1'
   goto, fine
endif
first=strpos(valore_all,'=')
last=strpos(valore_all,'/')
if last eq -1 then last=10000.
len=last-first
valore=strmid(valore_all,first+1,len-1)

fine:
end
