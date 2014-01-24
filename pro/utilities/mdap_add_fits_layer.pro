pro mdap_add_fits_layer,fitsfile,data,exten_no,exten_name,string_to_add

writefits,fitsfile,data,/append
h1=HEADFITS(fitsfile,EXTEN = exten_no)                                         
sxaddpar,h1,exten_name,string_to_add
modfits,fitsfile,0,h1,exten_no=exten_no

end
