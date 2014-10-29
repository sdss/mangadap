pro make_all

dir = '../high_level/results_rss/'
gal_name='ma001_145.125955+21.253809'
do_maps,dir,gal_name,2.88 
;stop
do_maps,dir,'ma002_164.44418+36.282683',3.35   
do_maps,dir,'ma002_208.04816+13.999926',7.32    
do_maps,dir,'ma002_ugc5124',9.27            
do_maps,dir,'ma003_163.98025+36.861503' ,10.25 
do_maps,dir,'ma003_207.87863+14.092194', 16.95   
do_maps,dir,'ma003_ngc2916',23.67               
do_maps,dir,'ma004_165.05043+36.387252',2.48   
do_maps,dir,'ma004_207.58111+14.140749',1.12   
do_maps,dir,'ma005_142.790030+22.746507',3.95    
do_maps,dir,'ma005_162.49460+36.415032',4.55    
do_maps,dir,'ma005_207.29175+13.347008',2.37    
do_maps,dir,'ma006_164.02366+36.960032',4.77     
do_maps,dir,'ma006_206.70605+14.400476',2.82    
do_maps,dir,'ma007_142.787803+20.916828',3.16    
do_maps,dir,'ma008_163.24608+37.613401',13.86    
do_maps,dir,'ma008_209.23090+14.142252',9.61
do_maps,dir,'ma008_cgcg122-022',6.77 


;plot_spectra,gal_name
end 

pro do_maps,dir,gal_name,Reff
;        !PATH
;restore, dir+gal_name+'/'+gal_name+'_mdap_session.idl'

;input = dir+gal_name+'/'+gal_name+'_high_level.fits'
input = dir+gal_name+'_high_level.fits'
res = file_test(input)
spatial_binning_str=mrdfits(input,4)
if res ne 1 then goto, skippa
;stop

;image = readfits(input,header,exten_no=0)
readfits_2d,input,image,x2d=x2d,y2d=y2d

binning=2
mdap_recovery_2d,input,binning,'vel',star_vel;,fits_file='vel1.fits'
mdap_recovery_2d,input,binning,'disp',star_sig;,fits_file='sig1.fits'
data=mrdfits(input,5)
xbin_str=data.x
ybin_str=data.y


binning=3
mdap_recovery_2d,input,binning,'vel',gas_vel;,fits_file='vel2.fits'
mdap_recovery_2d,input,binning,'OIII_5007_FLUX',gas_ha;,fits_file='OIII.fits'
mdap_recovery_2d,input,binning,'Ha_6563_FLUX',gas_OIII;,fits_file='Ha.fits'

binning=1
mdap_recovery_2d,input,binning,'hb',ind_hb;,fits_file='Hbeta.fits'
mdap_recovery_2d,input,binning,'mgb',ind_mgb;,fits_file='Mgb.fits'
mdap_recovery_2d,input,binning,'fe5270',ind_fe1
mdap_recovery_2d,input,binning,'fe5335',ind_fe2
ind_fe=0.5*(ind_fe1+ind_fe2)

;restore, dir+gal_name+'/'+gal_name+'_mdap_session.idl'


dist_bins_str = sqrt(xbin_str^2+ybin_str^2)
central_bin_str = where(dist_bins_str eq min(dist_bins_str))
central_bin_str = central_bin_str[0]

junk = where(dist_bins_str ge max(dist_bins_str) /2. )
v = min(dist_bins_str[junk])
mid_bin_str = where(abs(dist_bins_str - v[0]) le 0.001)
mid_bin_str=mid_bin_str[0]


outer_bin_str = where(dist_bins_str eq max(dist_bins_str))
outer_bin_str =outer_bin_str[0]
;stop

; !PATH = Expand_Path('+/usr/local/exelis/idl/')+ ':' + $
;         Expand_Path('+/home/UNI/coccatol/data/astron/') ;+ ':' + $
; !DIR  = '/usr/local/exelis/idl82'

set_plot, 'ps'
device, filename=gal_name+'.ps',/color,xsize=20,ysize=20,xoffset=.0,yoffset=4,bits_per_pixel=100

sx=.24
x0=0.1
y0=0.1
space=.055
top=0.65

mdap_place_image,2.5*alog10(image),x2d,y2d,position=[x0,y0+top,x0+sx,y0+sx+top],extra_inputs=['noerase=0','charsize=.8'],itable=54;,min=,max=max(2.5*alog10(image))
basic_colors, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   


indici = where(spatial_binning_str eq central_bin_str)
oplot,x2d[indici],y2d[indici],psym=6,color=red,symsize=.6
plots,xbin_str[central_bin_str],ybin_str[central_bin_str],psym=7,symsize=.8,thick=2,color=white

indici = where(spatial_binning_str eq mid_bin_str)
oplot,x2d[indici],y2d[indici],psym=6,color=green,symsize=.6
plots,xbin_str[mid_bin_str],ybin_str[mid_bin_str],psym=7,symsize=.8,thick=2
;stop
indici = where(spatial_binning_str eq outer_bin_str)
oplot,x2d[indici],y2d[indici],psym=6,color=yellow,symsize=.6
plots,xbin_str[outer_bin_str],ybin_str[outer_bin_str],psym=7,symsize=.8,thick=2

;stop

md1 = max(star_vel(where(finite(star_vel) eq 1)))
md2 = min(star_vel(where(finite(star_vel) eq 1)))
mm = median(star_vel)
md = max([abs(mm-md1) , abs(mm-md2)])*1.2
mdap_place_image,star_vel,x2d,y2d,position=[x0+space+sx,y0+top,x0+sx+space+sx,y0+sx+top],extra_inputs=['noerase=1','charsize=.8'],min=mm-md,max=mm+md,itable='mdap_color.tab'
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[mm-md,mm+md],ctable='mdap_color.tab',label='V stars [km/sec]',csize=.8;,/left

md1 = max(star_sig(where(finite(star_sig) eq 1)))+5
mdap_place_image,star_sig,x2d,y2d,position=[x0+sx+space+sx+space,y0+top,x0+sx+space+sx+space+sx,y0+sx+top],extra_inputs=['noerase=1','charsize=.8'],min=0,max=md1,itable='mdap_color.tab'
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[0,md1],ctable='mdap_color.tab',label='!4r!3 stars [km/sec]',csize=.8;,/left


md1 = max(gas_vel(where(finite(gas_vel) eq 1)))
md2 = min(gas_vel(where(finite(gas_vel) eq 1)))
mm = median(gas_vel)
md = max([abs(mm-md1) , abs(mm-md2)])*1.2
mdap_place_image,gas_vel,x2d,y2d,position=[x0,y0+top-2.*space-sx,x0+sx,y0+sx+top-2.*space-sx],extra_inputs=['noerase=1','charsize=.8'],min=mm-md,max=mm+md,itable='mdap_color.tab'
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[mm-md,mm+md],ctable='mdap_color.tab',label='V gas [km/sec]',csize=.8;,/left

ttt=gas_ha+gas_OIII
ttt = ttt[where(finite(ttt) eq 1 )]
ttt = 2.5*alog10(ttt)
ttt = ttt[where(finite(ttt) eq 1 )]
mdap_place_image,2.5*alog10(gas_ha),x2d,y2d,position=[x0+space+sx,y0+top-2.*space-sx,x0+sx+space+sx,y0+sx+top-2.*space-sx],extra_inputs=['noerase=1','charsize=.8'],itable=54,min=min(ttt),max=max(ttt);,max=max(2.5*alog10(gas_ha))
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[0,1],ctable=54,label='Halpha',csize=.8

mdap_place_image,2.5*alog10(gas_OIII),x2d,y2d,position=[x0+sx+space+sx+space,y0+top-2.*space-sx,x0+sx+space+sx+space+sx,y0+sx+top-2.*space-sx],extra_inputs=['noerase=1','charsize=.8'],itable=54,min=min(ttt),max=max(ttt);,max=max(2.5*alog10(gas_ha))
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[0,1],ctable=54,label='OIII 5007',csize=.8


md1 = max(ind_hb(where(finite(ind_hb) eq 1)))
md2 = min(ind_hb(where(finite(ind_hb) eq 1)))
mm = median(ind_hb)
md = max([abs(mm-md1) , abs(mm-md2)])*1.2
mdap_place_image,ind_hb,x2d,y2d,position=[x0,y0+top-4.*space-2.*sx,x0+sx,y0+sx+top-4.*space-2.*sx],extra_inputs=['noerase=1','charsize=.8'],min=mm-md,max=mm+md,itable='mdap_color.tab'
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[mm-md,mm+md],ctable='mdap_color.tab',label='Hbeta [ang]',csize=.8;,/left

md1 = max(ind_mgb(where(finite(ind_mgb) eq 1)))
md2 = min(ind_mgb(where(finite(ind_mgb) eq 1)))
mm = median(ind_mgb)
md = max([abs(mm-md1) , abs(mm-md2)])*1.2
mdap_place_image,ind_mgb,x2d,y2d,position=[x0+space+sx,y0+top-4.*space-2.*sx,x0+sx+space+sx,y0+sx+top-4.*space-2.*sx],extra_inputs=['noerase=1','charsize=.8'],min=mm-md,max=mm+md,itable='mdap_color.tab'
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[mm-md,mm+md],ctable='mdap_color.tab',label='Mgb [ang]',csize=.8;,/left

md1 = max(ind_fe(where(finite(ind_fe) eq 1)))
md2 = min(ind_fe(where(finite(ind_fe) eq 1)))
mm = median(ind_fe)
md = max([abs(mm-md1) , abs(mm-md2)])*1.2
mdap_place_image,ind_fe,x2d,y2d,position=[x0+sx+space+sx+space,y0+top-4.*space-2.*sx,x0+sx+space+sx+space+sx,y0+sx+top-4.*space-2.*sx],extra_inputs=['noerase=1','charsize=.8'],min=mm-md,max=mm+md,itable='mdap_color.tab'
place_cbar,min(x2d),min(y2d)-0.2*(max(y2d)-min(y2d)),max(x2d)-min(x2d),1.,[mm-md,mm+md],ctable='mdap_color.tab',label='<Fe> [ang]',csize=.8;,/left


;[x0+sx+space+sx+space,y0+.6-2.*space-sx,x0+sx+space+sx+space+sx,y0+sx+.6-2.*space-sx]
device, /close
set_plot, 'x'
print,gal_name +' done' 

skippa:
end
pro junk
x0=0.06
x1=0.98
y0=0.1
y1=0.98
ny=3
dy=(y1-y0)/ny
pos1=[x0,y0+0.*dy-0.03,x1,y0+1.*dy-0.03]
pos2=[x0,y0+1.*dy-0.03,x1,y0+2.*dy-0.03]
pos3=[x0,y0+2.*dy-0.03,x1,y0+3.*dy-0.03]
xrange=[3500,7300]

;window,0,retain=2,xsiz=1000,ysize=1200
;stop
indxx = where(WAVELENGTH_OUTPUT_str ge 4000 and WAVELENGTH_OUTPUT_str le 7000)
lrange=max(WAVELENGTH_OUTPUT_str[indxx])-min(WAVELENGTH_OUTPUT_str[indxx])
step = lrange/n_elements(indxx)

basic_colors, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   

plot, exp(log_wav_str),log_spc_str[*,central_bin_str],xrange=xrange,position=pos1,xtitle='angstrom',/norm,xstyle=1,yrange=[0,3.*median(log_spc_str[*,central_bin_str],/even)]
;oplot,exp(log_wav_str),log_err_str[*,central_bin_str],color=yellow
xyouts,0.1,0.34 ,'Dist (arcs): '+round_str(dist_bins_str[central_bin_str],2)+'; R/Re= '+round_str(dist_bins_str[central_bin_str]/reff,2),/norm
;xyouts,0.1,0.22 ,'S/N (per AA): '+round_str(bin_sn_str[central_bin_str],2),/norm
indici = where(spatial_binning_str eq central_bin_str)
xyouts,0.1,0.30 ,'Bin area (arcs!U2!N): '+round_str(n_elements(indici)*cdelt1*cdelt2,2),/norm
s = median(BEST_FIT_MODEL_STR[central_bin_str, indxx],/even)
n = robust_sigma(RESIDUALS_STR[central_bin_str, indxx])
xyouts,0.1,0.32 ,'S/N (per AA): '+round_str(s/n*sqrt(step),2),/norm

plot, exp(log_wav_str),log_spc_str[*,mid_bin_str],xrange=xrange,position=pos2,/norm,xtickname=replicate(' ',20),/noerase,xstyle=1,yrange=[0,3.*median(log_spc_str[*,mid_bin_str],/even)]
;oplot,exp(log_wav_str),log_err_str[*,mid_bin_str],color=yellow
xyouts,0.1,0.34+dy ,'Dist (arcs): '+round_str(dist_bins_str[mid_bin_str],2)+'; R/Re= '+round_str(dist_bins_str[mid_bin_str]/reff,2),/norm
;xyouts,0.1,0.22+dy ,'S/N (per AA): '+round_str(bin_sn_str[mid_bin_str],2),/norm
indici = where(spatial_binning_str eq mid_bin_str)
xyouts,0.1,0.30+dy ,'Bin area (arcs!U2!N): '+round_str(n_elements(indici)*cdelt1*cdelt2,2),/norm
s = median(BEST_FIT_MODEL_STR[mid_bin_str, indxx],/even)
n = robust_sigma(RESIDUALS_STR[mid_bin_str, indxx])
xyouts,0.1,0.32+dy ,'S/N (per AA): '+round_str(s/n*sqrt(step),2),/norm

plot, exp(log_wav_str),log_spc_str[*,outer_bin_str],xrange=xrange,position=pos3,/norm,xtickname=replicate(' ',20),/noerase,xstyle=1,title=gal_name,yrange=[0,3.*median(log_spc_str[*,outer_bin_str],/even)]
;oplot,exp(log_wav_str),log_err_str[*,outer_bin_str],color=yellow
xyouts,0.1,0.34+2.*dy ,'Dist (arcs): '+round_str(dist_bins_str[outer_bin_str],2)+'; R/Re= '+round_str(dist_bins_str[outer_bin_str]/reff,2),/norm
;xyouts,0.1,0.22+2*dy ,'S/N (per AA): '+round_str(bin_sn_str[outer_bin_str],2),/norm
indici = where(spatial_binning_str eq outer_bin_str)
xyouts,0.1,0.30+2*dy ,'Bin area (arcs!U2!N): '+round_str(n_elements(indici)*cdelt1*cdelt2,2),/norm
s = median(BEST_FIT_MODEL_STR[outer_bin_str, indxx],/even)
n = robust_sigma(RESIDUALS_STR[outer_bin_str, indxx])
xyouts,0.1,0.32+2.*dy ,'S/N (per AA): '+round_str(s/n*sqrt(step),2),/norm










basic_colors, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   

plot, exp(log_wav_str),log_spc_str[*,central_bin_str],xrange=xrange,position=pos1,xtitle='angstrom',/norm,xstyle=1,yrange=[0,3.*median(log_spc_str[*,central_bin_str],/even)]
oplot,wavelength_output_str,best_fit_model_str[central_bin_str,*],color=red
xyouts,0.1,0.34 ,'Dist (arcs): '+round_str(dist_bins_str[central_bin_str],2)+'; R/Re= '+round_str(dist_bins_str[central_bin_str]/reff,2),/norm
;xyouts,0.1,0.22 ,'S/N (per AA): '+round_str(bin_sn_str[central_bin_str],2),/norm
indici = where(spatial_binning_str eq central_bin_str)
xyouts,0.1,0.30 ,'Bin area (arcs!U2!N): '+round_str(n_elements(indici)*cdelt1*cdelt2,2),/norm
s = median(BEST_FIT_MODEL_STR[central_bin_str, indxx],/even)
n = robust_sigma(RESIDUALS_STR[central_bin_str, indxx])
xyouts,0.1,0.32 ,'S/N (per AA): '+round_str(s/n*sqrt(step),2),/norm

plot, exp(log_wav_str),log_spc_str[*,mid_bin_str],xrange=xrange,position=pos2,/norm,xtickname=replicate(' ',20),/noerase,xstyle=1,yrange=[0,3.*median(log_spc_str[*,mid_bin_str],/even)]
oplot,wavelength_output_str,best_fit_model_str[mid_bin_str,*],color=red
;oplot,exp(log_wav_str),log_err_str[*,mid_bin_str],color=yellow
xyouts,0.1,0.34+dy ,'Dist (arcs): '+round_str(dist_bins_str[mid_bin_str],2)+'; R/Re= '+round_str(dist_bins_str[mid_bin_str]/reff,2),/norm
;xyouts,0.1,0.22+dy ,'S/N (per AA): '+round_str(bin_sn_str[mid_bin_str],2),/norm
indici = where(spatial_binning_str eq mid_bin_str)
xyouts,0.1,0.30+dy ,'Bin area (arcs!U2!N): '+round_str(n_elements(indici)*cdelt1*cdelt2,2),/norm
s = median(BEST_FIT_MODEL_STR[mid_bin_str, indxx],/even)
n = robust_sigma(RESIDUALS_STR[mid_bin_str, indxx])
xyouts,0.1,0.32+dy ,'S/N (per AA): '+round_str(s/n*sqrt(step),2),/norm

plot, exp(log_wav_str),log_spc_str[*,outer_bin_str],xrange=xrange,position=pos3,/norm,xtickname=replicate(' ',20),/noerase,xstyle=1,title=gal_name,yrange=[0,3.*median(log_spc_str[*,outer_bin_str],/even)]
oplot,wavelength_output_str,best_fit_model_str[outer_bin_str,*],color=red
;oplot,exp(log_wav_str),log_err_str[*,outer_bin_str],color=yellow
xyouts,0.1,0.34+2.*dy ,'Dist (arcs): '+round_str(dist_bins_str[outer_bin_str],2)+'; R/Re= '+round_str(dist_bins_str[outer_bin_str]/reff,2),/norm
;xyouts,0.1,0.22+2*dy ,'S/N (per AA): '+round_str(bin_sn_str[outer_bin_str],2),/norm
indici = where(spatial_binning_str eq outer_bin_str)
xyouts,0.1,0.30+2*dy ,'Bin area (arcs!U2!N): '+round_str(n_elements(indici)*cdelt1*cdelt2,2),/norm
s = median(BEST_FIT_MODEL_STR[outer_bin_str, indxx],/even)
n = robust_sigma(RESIDUALS_STR[outer_bin_str, indxx])
xyouts,0.1,0.32+2.*dy ,'S/N (per AA): '+round_str(s/n*sqrt(step),2),/norm








save,filename=gal_name+'session.idl',image,star_vel,x2d,y2d,log_spc_str,log_wav_str,best_fit_model_str,wavelength_output;,/variables
device, /close
set_plot, 'x'
print,gal_name +' done' 

skippa:
end
