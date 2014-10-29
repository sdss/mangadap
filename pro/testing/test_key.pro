
FUNCTION check_data, data=data
;	if keyword_set(data) then return, 1
	if n_elements(data) then return, 1
	return, 0
END


PRO set_data, data=data
    if check_data(data=data) eq 1 then begin
	data=10.0
    endif else $
	data=0.0
    print, data

END



