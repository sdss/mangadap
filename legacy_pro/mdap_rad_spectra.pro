pro make_all

dir = '../results_datacubes/'
gal_name='ma001_145.125955+21.253809'
do_rad_spectra,dir,gal_name,2.88 
;stop
do_rad_spectra,dir,'ma002_164.44418+36.282683',3.35   
do_rad_spectra,dir,'ma002_208.04816+13.999926',7.32    
do_rad_spectra,dir,'ma002_ugc5124',9.27            
do_rad_spectra,dir,'ma003_163.98025+36.861503' ,10.25 
do_rad_spectra,dir,'ma003_207.87863+14.092194', 16.95   
do_rad_spectra,dir,'ma003_ngc2916',23.67               
do_rad_spectra,dir,'ma004_165.05043+36.387252',2.48   
do_rad_spectra,dir,'ma004_207.58111+14.140749',1.12   
do_rad_spectra,dir,'ma005_142.790030+22.746507',3.95    
do_rad_spectra,dir,'ma005_162.49460+36.415032',4.55    
do_rad_spectra,dir,'ma005_207.29175+13.347008',2.37    
do_rad_spectra,dir,'ma006_164.02366+36.960032',4.77     
do_rad_spectra,dir,'ma006_206.70605+14.400476',2.82    
do_rad_spectra,dir,'ma007_142.787803+20.916828',3.16    
do_rad_spectra,dir,'ma008_163.24608+37.613401',13.86    
do_rad_spectra,dir,'ma008_209.23090+14.142252',9.61
do_rad_spectra,dir,'ma008_cgcg122-022',6.77 


end
pro do_rad_spectra,dir,gal_name,Reff

input = dir+gal_name+'/'+gal_name+'_mdap_session.idl'
res = file_test(input)
if res ne 1 then goto, skippa

restore,input

xrange=[3600,10200]
ny=6
x0=0.1
x1=0.98
y0=0.1
y1=0.98
dy=(y1-y0)/ny
pos1=[x0,y0+0.0*dy,x1,y0+1.0*dy]
pos2=[x0,y0+1.0*dy,x1,y0+2.0*dy]
pos3=[x0,y0+2.0*dy,x1,y0+3.0*dy]
pos4=[x0,y0+3.0*dy,x1,y0+4.0*dy]
pos5=[x0,y0+4.0*dy,x1,y0+5.0*dy]
pos6=[x0,y0+5.0*dy,x1,y0+6.0*dy]
set_plot, 'ps'
device, filename=gal_name+'_spectra.ps',xsize=15,ysize=25,yoffset=1,xoffset=3

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,0],xrange=xrange,$
     position=pos1,xstyle=1,xtitle='wavelength [angstrom]'
xyouts,x0+0.05,0.9/6.*1+.06,'BIN 1:'+mdap_stc(r_bin[0]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,1],xrange=xrange,$
     position=pos2,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*2+.06,'BIN 2:'+mdap_stc(r_bin[1]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,2],xrange=xrange,$
     position=pos3,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*3+.06,'BIN 3:'+mdap_stc(r_bin[2]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,3],xrange=xrange,$
     position=pos4,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*4+.06,'BIN 4:'+mdap_stc(r_bin[3]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,4],xrange=xrange,$
     position=pos5,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*5+.06,'BIN 5:'+mdap_stc(r_bin[4]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,5],xrange=xrange,$
     position=pos6,xstyle=1,xtickname=replicate(' ',20),title=gal_name,/noerase
xyouts,x0+0.05,0.9/6.*6+.06,'BIN 6:'+mdap_stc(r_bin[5]),/norm


xrange=[8000,10200]

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,0],xrange=xrange,$
     position=pos1,xstyle=1,xtitle='wavelength [angstrom]'
xyouts,x0+0.05,0.9/6.*1+.06,'BIN 1:'+mdap_stc(r_bin[0]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,1],xrange=xrange,$
     position=pos2,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*2+.06,'BIN 2:'+mdap_stc(r_bin[1]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,2],xrange=xrange,$
     position=pos3,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*3+.06,'BIN 3:'+mdap_stc(r_bin[2]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,3],xrange=xrange,$
     position=pos4,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*4+.06,'BIN 4:'+mdap_stc(r_bin[3]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,4],xrange=xrange,$
     position=pos5,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*5+.06,'BIN 5:'+mdap_stc(r_bin[4]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,5],xrange=xrange,$
     position=pos6,xstyle=1,xtickname=replicate(' ',20),title=gal_name,/noerase
xyouts,x0+0.05,0.9/6.*6+.06,'BIN 6:'+mdap_stc(r_bin[5]),/norm


xrange=[3600, 7000]

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,0],xrange=xrange,$
     position=pos1,xstyle=1,xtitle='wavelength [angstrom]'
xyouts,x0+0.05,0.9/6.*1+.06,'BIN 1:'+mdap_stc(r_bin[0]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,1],xrange=xrange,$
     position=pos2,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*2+.06,'BIN 2:'+mdap_stc(r_bin[1]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,2],xrange=xrange,$
     position=pos3,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*3+.06,'BIN 3:'+mdap_stc(r_bin[2]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,3],xrange=xrange,$
     position=pos4,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*4+.06,'BIN 4:'+mdap_stc(r_bin[3]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,4],xrange=xrange,$
     position=pos5,xstyle=1,xtickname=replicate(' ',20),/noerase
xyouts,x0+0.05,0.9/6.*5+.06,'BIN 5:'+mdap_stc(r_bin[4]),/norm

plot, exp(loglam_gal_rbin), radially_binned_spectra[*,5],xrange=xrange,$
     position=pos6,xstyle=1,xtickname=replicate(' ',20),title=gal_name,/noerase
xyouts,x0+0.05,0.9/6.*6+.06,'BIN 6:'+mdap_stc(r_bin[5]),/norm


device, /close
set_plot, 'x'

skippa:
end
