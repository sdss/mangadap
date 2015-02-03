pro mdap_get_2d_map,values,binning,map

   values_ = reform(values)
   map = binning*0./0.
   for i = 0, max(binning) do begin
      ind = where(binning eq i)
      if ind[0] ne -1 then map[ind] = values_[i]
   endfor



end
