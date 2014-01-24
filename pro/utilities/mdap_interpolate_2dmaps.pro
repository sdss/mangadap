pro mdap_interpolate_2dmaps,input,x,y,x2d_full,y2d_full,x_out,y_out,output

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

sz=size(result)
x_out_=interpol(findgen(sz[1]),range(startx,endx,nx),x_out)
y_out_=interpol(findgen(sz[2]),range(starty,endy,ny),y_out)
for i = 0, n_elements(OUTPUT)-1 do output[i] = bilinear(result,x_out_[i],y_out_[i])

;triangulate,x,y,tr,b
;output = TRIGRID( X, Y, input, Tr)
;stop
end
