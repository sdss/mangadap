#!/bin/bash

if [ "$#" -ne 1 ]; then
    dir=$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
else
    dir=${1}
fi
echo ${dir}

((tot=0))
((fin=0))
for donef in $( ls */*/*.done ); do
    root=$( echo $donef | cut -f 1 -d '.')
    ((++tot))
    if [ $( grep -c SUCCESSFULLY ${root}.out ) -eq 1 ]; then 
        echo ${root}
        ((++fin))
    fi
done

echo "Successfully processed ${fin} of ${tot} files."

exit


