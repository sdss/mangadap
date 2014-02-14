pro mdap_interpolate_2dmaps,input,x,y,x2d_full,y2d_full,x_out,y_out,output,full_grid=full_grid

;
; PURPOSE
; This routine: 
;    i) uses the GRID_TPS IDL function to interpolate a set of N input
;       values, defined over an irregular grid, on to a regular N'xM' grid.
;   ii) interpolates the values interpolated over the regular grid
;       (point i) and re-interpolates them over an irregular grid using
;       the IDL BILINEAR function.
;
;
;INPUTS
;
; input   N elements array containing the values to be interpolated. 
;
; x       N elements array containing the X coordinates where the N input values are defined.
;
; y       N elements array containing the Y coordinates where the N input values are defined.
;
; x2d_full N'xM' elements array containing the X coordiantes of a regularly spaced grid.
;
; y2d_full N'xM' elements array containing the Y coordiantes of a regularly spaced grid.
;
; x_out M elements array containing the X coordinates where the input values need to be interpolated
;
; y_out M elements array containing the X coordinates where the input values need to be interpolated
;
;
; OUTPUTS
;
; output M elements array containing the interpolated input values on the x_out,y_out coordiantes.
;
;
; OPTIONAL OUTPUT
;
; full_grid N'xM' elements array containing the interpolated input values on the regular grid with coordinates x2d_full, y2d_full 
;
;


index_good = where(finite(input) eq 1)
if n_elements(index_good) le 7 then index_good = indgen(n_elements(x))

sz=size(x2d_full)
output = dblarr(n_elements(x_out),/nozero)
;output = min_curve_surf(input,x,y,xout=x_out,yout=y_out)

startx = min(x2d_full)-1
starty = min(y2d_full)-1
endx = max(x2d_full)+1
endy = max(y2d_full)+1
nx = sz[1]*2.+1
ny = sz[2]*2.+1
dx = (endx - startx)/(nx-1) 
dy = (endy - starty)/(ny-1) 
result = GRID_TPS (X[index_good], Y[index_good], input[index_good],ngrid=[nx,ny],start=[startx,starty],delta=[dx,dy])
full_grid=result
sz=size(result)
x_out_=interpol(findgen(sz[1]),range(startx,endx,nx),x_out)
y_out_=interpol(findgen(sz[2]),range(starty,endy,ny),y_out)
for i = 0, n_elements(OUTPUT)-1 do output[i] = bilinear(result,x_out_[i],y_out_[i])

;triangulate,x,y,tr,b
;output = TRIGRID( X, Y, input, Tr)
;stop
end
