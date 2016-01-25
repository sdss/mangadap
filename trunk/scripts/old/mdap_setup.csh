#!/bin/csh -f

# Check for proper usage
if ($#argv != 3) then
	echo "ERROR: Incorrect number of arguments"
	echo "USAGE: source setup.sh <output file> <type? (BOTH, CUBE, or RSS)> <overwrite? (1-yes;0-no)>"
	exit 1
endif

# File exists, but overwrite not allowed
if ( -e $1 && $3 == 0) then
	echo "$1 already exists! Set overwrite flag or remove."
	echo "USAGE: source setup.sh <output file> <type? (BOTH, CUBE, or RSS)> <overwrite? (1-yes;0-no)>"
	exit 0
endif

# Set the environmental variables that will be read by IDL
source {$MDAP_DEV}/scripts/mdap_environment.sh

# Check for the existence of the analysis directories and create them if
# necessary
# if ( ! -e ${MANGA_SPECTRO}/analysis ) then
# 	mkdir ${MANGA_SPECTRO}/analysis
# endif

if ( ! -e ${MANGA_SPECTRO}/analysis/${VERSDAP} ) then
	mkdir -p ${MANGA_SPECTRO}/analysis/${VERSDAP}
endif

# Find files in redux/*/stack
rm -f $1
if ( $2 == 'BOTH') then
    set list='*/stack/*.fits*'
else if ( $2 == 'CUBE' ) then
    set list='*/stack/*CUBE.fits*'
else if ( $2 == 'RSS' ) then
    set list='*/stack/*RSS.fits*'
else
    echo "ERROR: Unknown file type: $2. Should be (BOTH, CUBE, or RSS)"
    exit 1
endif

printf '' > $1
foreach file ( `ls {$MANGA_SPECTRO_REDUX}/$list` )
    set ifun = ` echo $file | cut -f 3 -d '-' `
    if ( $ifun > 1900 ) then
        echo $file >> $1
    endif
end

# TODO: Check for associated files produced by the DAP

# Report
set nf = ` wc $1 | awk '{ print $1 }' `
echo "Number of files found to process: $nf"

exit 0



