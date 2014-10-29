pro make_all

dir = '../high_level/results_datacubes/'
gal_name='ma001_145.125955+21.253809'
do_ind_profiles,dir,gal_name,2.88 
;stop
do_ind_profiles,dir,'ma002_164.44418+36.282683',3.35   
do_ind_profiles,dir,'ma002_208.04816+13.999926',7.32    
do_ind_profiles,dir,'ma002_ugc5124',9.27            
do_ind_profiles,dir,'ma003_163.98025+36.861503' ,10.25 
do_ind_profiles,dir,'ma003_207.87863+14.092194', 16.95   
do_ind_profiles,dir,'ma003_ngc2916',23.67               
do_ind_profiles,dir,'ma004_165.05043+36.387252',2.48   
do_ind_profiles,dir,'ma004_207.58111+14.140749',1.12   
do_ind_profiles,dir,'ma005_142.790030+22.746507',3.95    
do_ind_profiles,dir,'ma005_162.49460+36.415032',4.55    
do_ind_profiles,dir,'ma005_207.29175+13.347008',2.37    
do_ind_profiles,dir,'ma006_164.02366+36.960032',4.77     
do_ind_profiles,dir,'ma006_206.70605+14.400476',2.82    
do_ind_profiles,dir,'ma007_142.787803+20.916828',3.16    
do_ind_profiles,dir,'ma008_163.24608+37.613401',13.86    
do_ind_profiles,dir,'ma008_209.23090+14.142252',9.61
do_ind_profiles,dir,'ma008_cgcg122-022',6.77 


;plot_spectra,gal_name
end 
pro do_ind_profiles,dir,gal_name,reff
;window,0,retain=2,xsize=800,ysize=1200
input = dir+gal_name+'_high_level.fits'
res = file_test(input)
if res ne 1 then goto, skippa

set_plot, 'ps'
device, filename=gal_name+'_ind_rad_prof.ps',xsize=10,ysize=20,xoffset=2,yoffset=1

data=mrdfits(input,9,/silent)
data.amaj=data.amaj/Reff
ny = 8
x0=0.12
x1=0.98
y0=0.1
y1=0.97
dy = (y1-y0)/ny

pos1=[x0,y0+0.0*dy,x1,y0+1.0*dy]
pos2=[x0,y0+1.0*dy,x1,y0+2.0*dy]
pos3=[x0,y0+2.0*dy,x1,y0+3.0*dy]
pos4=[x0,y0+3.0*dy,x1,y0+4.0*dy]
pos5=[x0,y0+4.0*dy,x1,y0+5.0*dy]
pos6=[x0,y0+5.0*dy,x1,y0+6.0*dy]
pos7=[x0,y0+6.0*dy,x1,y0+7.0*dy]
pos8=[x0,y0+7.0*dy,x1,y0+8.0*dy]

xrange=[0,max(data.amaj_up/Reff)+Reff/20.]
y=data.D4000_rad
ey=data.D4000_rad_err

yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos1,ytitle='D4000',xtitle='a [arcsec/Reff]',xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey

y=data.Hb_rad
ey=data.Hb_rad_err
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos2,xtickname=replicate(' ',20),ytitle='Hbeta',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey


y=data.Mgb_rad
ey=data.Mgb_rad_err
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos3,xtickname=replicate(' ',20),ytitle='Mgb',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey

y=data.Fe5270_rad
ey=data.Fe5270_rad_err
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos4,xtickname=replicate(' ',20),ytitle='Fe5270',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey

y=data.Fe5335_rad
ey=data.Fe5335_rad_err
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos5,xtickname=replicate(' ',20),ytitle='Fe5335',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey

y=data.NaD_rad
ey=data.NaD_rad_err
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos6,xtickname=replicate(' ',20),ytitle='NaD',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey




y=(data.CaII0p86a_rad+data.CaII0p86b_rad+data.CaII0p86c_rad)/3.
ey=sqrt(data.CaII0p86a_rad_err^2+data.CaII0p86b_rad_err^2+data.CaII0p86c_rad_err^2)/sqrt(3)
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos7,xtickname=replicate(' ',20),ytitle='CaII0.86',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1
oploterror,data.amaj,y,data.amaj*0.,ey


y=data.NaI0p82_rad
ey=data.NaI0p82_rad_err
yrange=[min(y-ey)-0.05,max(y+ey)+0.05]
plot,data.amaj,y,psym=4,pos=pos8,xtickname=replicate(' ',20),ytitle='NaI0.82',/noerase,xrange=xrange,xstyle=1,$
yrange=yrange,ystyle=1,title=gal_name
oploterror,data.amaj,y,data.amaj*0.,ey



device, /close
set_plot, 'x'
skippa:
end
