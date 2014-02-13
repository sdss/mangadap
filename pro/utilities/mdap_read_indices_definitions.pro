pro mdap_read_indices_definitions,absorption_line_indices,indices=indices

readcol,absorption_line_indices,id,name,passband0,passband1,blue_cont0,blue_cont1,red_cont0,red_cont1,units,$
                                                format='I,A,D,D,D,D,D,D,A',comment='#',/silent                                                        

;; not correct: 
;passband_TiO2=[6191.375,6273.875]                       
;blue_cont_TiO2=[6068.375,6143.375]                      
;red_cont_TiO2=[6374.375,6416.875]                       


indices=replicate({id:0, name:'n',units:'ang',$
                  passband:[0.d,0.d], red_cont:[0.d,0.d],blue_cont:[0.d,0.d]}, n_elements(id))


for i = 0, n_elements(id)-1 do begin
   indices[i].id = id[i]
   indices[i].name = name[i]
   indices[i].units=units[i]
   indices[i].passband=[passband0[i],passband1[i]]
   indices[i].blue_cont=[blue_cont0[i],blue_cont1[i]]
   indices[i].red_cont=[red_cont0[i],red_cont1[i]]
endfor


end
