function mdap_round_str,valore,decimals
 str=mdap_stc(round(valore*(10.^decimals))/10.^decimals)
 pos=strpos(str,'.')
 res=strmid(str,0,pos+decimals+1)
return, res
end
