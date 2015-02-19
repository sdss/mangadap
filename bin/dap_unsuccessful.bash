#!/bin/bash

if [ "$#" -ne 1 ]; then
    dir=$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER
else
    dir=${1}
fi
echo ${dir}

((tot=0))
((flt=0))
for out in $( ls */*/*.out ); do
    ((++tot))
    if [ $( grep -c SUCCESSFULLY ${out} ) -ne 1 ]; then 
        echo ${out}
        ((++flt))
    fi
done

echo "Fault in ${flt} of ${tot} files."

exit


