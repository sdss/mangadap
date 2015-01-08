#!/bin/bash

# Check for proper usage
if [ $# -ne 3 ]; then
	echo "ERROR: Incorrect number of arguments"
	echo "USAGE: source mdap_setup.bash <output file> <type? (BOTH, CUBE, or RSS)> <overwrite? (1-yes;0-no)>"
	return 1
fi

# File exists, but overwrite not allowed
if [ -e $1 ] && [ $3 -eq 0 ]; then
	echo "$1 already exists! Set overwrite flag or remove."
	echo "USAGE: source mdap_setup.bash <output file> <type? (BOTH, CUBE, or RSS)> <overwrite? (1-yes;0-no)>"
	return 0
fi

# Set the environmental variables that will be read by IDL
source $MDAP_DEV/scripts/mdap_environment.bash

# Check for the existence of the analysis directories and create them if
if [ ! -e $MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER ]; then
	mkdir -p $MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER
fi

# Find files in redux/*/stack
rm -f $1
if [ $2 == 'BOTH' ]; then
    list='*/stack/*.fits*'
elif [ $2 == 'CUBE' ]; then
    list='*/stack/*CUBE.fits*'
elif [ $2 == 'RSS' ]; then
    list='*/stack/*RSS.fits*'
else
    echo "ERROR: Unknown file type: $2. Should be (BOTH, CUBE, or RSS)"
    return 1
fi

printf '' > $1
for file in $(ls $MANGA_SPECTRO_REDUX/$list) ; do
    ifun=( ` echo $file | cut -f 3 -d '-' ` )
    if [ $ifun -gt 1900 ]; then
        echo $file >> $1
    fi
done

# TODO: Check for associated files produced by the DAP

# Report
nf=( ` wc $1 | awk '{ print $1 }' ` )
echo "Number of files found to process: $nf"

return 0



