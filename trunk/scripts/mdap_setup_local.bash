#!/bin/bash

# Setup for running the MaNGA DAP using a local installation with local
# paths

# Ignores minibundles

# Check for proper usage
if [ $# -ne 4 ]; then
	echo "ERROR: Incorrect number of arguments"
	echo "USAGE: mdap_setup_local.bash <file path> <output file> <type? (BOTH, CUBE, or RSS)> <overwrite? (1-yes;0-no)>"
	exit 1
fi

# File exists, but overwrite not allowed
if [ -e $2 ] && [ $4 -eq 0 ]; then
	echo "$2 already exists! Set overwrite flag or remove."
	echo "USAGE: mdap_setup_local.bash <file path> <output file> <type? (BOTH, CUBE, or RSS)> <overwrite? (1-yes;0-no)>"
	exit 0
fi

# Find files in file_path
rm -f $2
if [ $3 == 'BOTH' ]; then
    list='*.fits*'
elif [ $3 == 'CUBE' ]; then
    list='*CUBE.fits*'
elif [ $3 == 'RSS' ]; then
    list='*RSS.fits*'
else
    echo "ERROR: Unknown file type: $3. Should be (BOTH, CUBE, or RSS)"
    exit 1
fi

printf '' > $2
for file in $(ls ${1}/$list) ; do
    ifun=( ` echo $file | cut -f 3 -d '-' ` )
    if [ $ifun -gt 1900 ]; then
        echo $file >> $2
    fi
done

# Report
nf=( ` wc $2 | awk '{ print $2 }' ` )
echo "Number of files found to process: $nf"

exit 0



