function mdap_stc,value,integer=integer

if not keyword_set(integer) then s = strcompress(string(value),/remove_all)
if keyword_set(integer) then s=strcompress(string(long(value)),/remove_all)

return,s
end
