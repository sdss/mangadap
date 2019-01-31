#!/bin/bash

if [ "$#" -ne 1 ]; then
    dir=$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
else
    dir=${1}
fi
echo ${dir}

((fin=0))
((run=0))
for strtf in $( ls */*/*.started ); do
    root=$( echo $strtf | cut -f 1 -d '.')
    if [ -e ${root}.done ]; then
        ((++fin))
    else
        echo ${root}
        ((++run))
    fi
done

echo "Completed ${fin}; ${run} currently running."

exit


