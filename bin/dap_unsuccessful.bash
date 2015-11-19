#!/bin/bash

if [ "$#" -ne 1 ]; then
    dir=$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
else
    dir=${1}
fi
echo ${dir}

((tot=0))
((flt=0))
for donef in $( ls ${dir}/*/*/*.done ); do
    root=$( echo $donef | cut -f 1 -d '.')
    ((++tot))
    if [ $( grep -c SUCCESSFULLY ${root}.out ) -ne 1 ]; then 
        echo ${root}
        ((++flt))
    fi
done

echo "Fault in ${flt} of ${tot} files."

exit


