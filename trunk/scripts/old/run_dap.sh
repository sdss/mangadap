#!/bin/bash
################################################################################
# Revision History:
#	13 Nov 2014: (KBW) Original implementation
#
# Description
#
#	Execute a single instance of the MANGA DAP using the provided setup
#	procedure and the root name of a designated DRP-produced fits file to
#	process.
#
#	The script:
#
#	    - does some basic checks of the directory structure and the
#	      availability of the DRP file.
#
#	    - gets the list of available DRP files, excluding the mini-bundles,
#	      hard-wired to be written to 
#	    
#	    list_file="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/tmp_run_dap.lst"
#
#	    - checks this file against the existing list of available files, if
#	      it exists; this file is hard-wired to be
#
#	    drp_list="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/drp_completed.lst"
#
#	    - if $list_file and $drp_list are different, or the required input
#	      table for the DAP, hard-wired to be
#
#	    input_table="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/manga_dap_table.inp"
#
#	      is not available, it is created using the IDL procedure
#	      MDAP_DRPALL_INFO.  Eventually this should gather the required
#	      information from the DRPall file.  For now, this procedure gathers
#	      the information directly from the NSA catalog by first matching
#	      the object coordinates in the DRP fits file to the galaxies in the
#	      NSA catalog (MDAP_MATCH_OBS_NSA) and then creating the input table
#	      by pulling the needed information from the NSA catalog
#	      (MDAP_CREATE_INPUT_TABLE).  The match results are hard-wired to be output to
#
#	    nsa_match="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/manga_drp_nsa_match.inp"
#
#	      and the NSA catalog is hard-wired to 
#
#	    nsa_cat="$MANGAWORK_DIR/manga/target/temp/12-nsa_v1b_0_0_v2.fits.gz"
#
#	    - copies the provided setup procedure (the first argument) to
#	      $MANGADAP_DIR/pro/usr/mdap_execution_setup.pro so that it is
#	      accessible via IDL_PATH
#
#	    - determines which line in the provided input table has the
#	      information for the requested file to process
#
#	    - executes the DAP for that line of the input table
#
# Parameters are:
#	$1 - the IDL procedure file that is required for the setup of the DAP
#		see $MANGADAP_DIR/examples/mdap_setup.pro as an example
#	     So that this procedure can be viewed by IDL_PATH, it is moved to 
#		$MANGADAP_DIR/pro/usr/mdap_execution_setup.pro
#
#	$2 - the DRP fits file to analyze WITHOUT the fits.gz extension.  The
#	     file MUST be found at
#		$MANGA_SPECTRO_REDUX/$MANGADRP_VER/[plate]/stack/${2}.fits.gz
#	     where [plate] is drawn directly from the file name, which is
#	     expected to be of the form
#
#		manga-[PLATE]-[IFUDESIGN]-LOGCUBE.fits.gz
#
#	     as per the TRM.
#
#	$3 - spool file for the IDL output
#
# Usage example: run_dap.sh mdap_setup.pro manga-7495-12703-LOGCUBE spool_dap
#
################################################################################


################################################################################
# FUNCTIONS
read_file_element () {
	local c
	c="{if(NR==$2)print\$$3}"
	echo "$(awk $c $1)"
}
read_file_column () {
	local n=0
	local c
	local d
	c="{if(NR>$2)print\$$3}"
	for entry in $(awk $c $1); do
	    d[$n]=$entry
	    ((n++))
	done
	echo "${d[@]}"
}
array_size () {
	local passed_array
	passed_array=( `echo "$1"` )
	echo ${#passed_array[@]}
}
# Take the difference $1 - $2
diff() {
	local out
	if [ $(echo "$2<=0" | bc -l) -eq 1 ]; then
	    out=( `echo "$1+sqrt($2*$2)" | bc -l` )
	else
	    out=( `echo "$1-$2" | bc -l` )
	fi
	echo $out
}
# Usage, e.g.: if [ $(find_string $1 "hos") -eq 1 ]; then
find_string () {
	case $1 in
	    # match exact string
	    "$2") echo 1 ;;
	    # match start of string
	    "$2"*) echo 1 ;;
	    # match end of string
	    *"$2") echo 1 ;;
	    # string $2 can be anywhere in string $1
	    *"$2"*) echo 1 ;;
	    # not matched
	    *) echo 0 ;;
	esac
}
################################################################################

# Start time
date

################################################################################
# Check the number of arguments
if [ $# -ne 3 ]; then
    echo "Usage: run_dap.sh <setup procedure file> <DRP fits file> <IDL spool file>"
    exit 1
fi

echo "Running script with parameters: $1 $2 $3"

################################################################################
# Check that the setup file exists
if [ ! -e ${1} ]; then
    echo "ERROR: File does not exist: ${1}"
    exit 1
fi

################################################################################
# Check that the DRP file exists
plate=( `echo ${2} | cut -d '-' -f 2` )
if [ ! -e $MANGA_SPECTRO_REDUX/$MANGADRP_VER/${plate}/stack/${2}.fits.gz ]; then
    echo "ERROR: File does not exist: $MANGA_SPECTRO_REDUX/$MANGADRP_VER/${plate}/stack/${2}.fits.gz"
    exit 1
fi

################################################################################
# Make sure that the main output directory exists
if [ ! -d $MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER ]; then
    echo "ERROR: Output directory does not exist: $MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER"
    exit 1
fi

################################################################################
# Update the list of available DRP files
type='BOTH'
if [ ${type} == 'BOTH' ]; then
    list='*/stack/*.fits*'
elif [ ${type} == 'CUBE' ]; then
    list='*/stack/*CUBE.fits*'
elif [ ${type} == 'RSS' ]; then
    list='*/stack/*RSS.fits*'
else
    echo "ERROR: Unknown file type: ${type}.  Should be BOTH, CUBE or RSS."
    exit 1
fi

list_file="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/tmp_run_dap.lst"
rm -f $list_file
#echo $list_file
printf '' > $list_file
for file in $(ls $MANGA_SPECTRO_REDUX/$MANGADRP_VER/$list) ; do
    ifun=( ` echo $file | cut -d '-' -f 3 ` )
    # Do not consider the mini-bundles
    if [ $ifun -gt 1900 ]; then
	echo $file >> $list_file
    fi
done

drp_list="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/drp_completed.lst"
#echo $drp_list
((build_input_table=0))
if [ ! -e $drp_list ]; then
    ((build_input_table=1))
else
    cmp -s $list_file $drp_list >/dev/null
    if [ $? -eq 1 ]; then
	((build_input_table=1))
    fi
fi

# The file list did not previously exist or has changed
input_table="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/manga_dap_table.inp"
if [ ${build_input_table} -eq 1 ] || [ ! -e $input_table ]; then
    echo "Building/Updating DAP input table ..."
    # nsa_match="$MANGA_SPECTRO_ANALYSIS/$MANGADAP_VER/manga_drp_nsa_match.inp"

    mv $list_file $drp_list
    cmd='idl -quiet -e mdap_drpall_info'
    echo `$cmd`
else
    rm -f $list_file
fi

# Move the setup file so that it's accessible by IDL_PATH
if [ -e $MANGADAP_DIR/pro/usr/mdap_execution_setup.pro ]; then
    echo "WARNING: Removing old mdap_execution_setup.pro file!"
    rm -f $MANGADAP_DIR/pro/usr/mdap_execution_setup.pro
fi
cp $1 $MANGADAP_DIR/pro/usr/mdap_execution_setup.pro

# Get the index of the requested file to analyze
declare -a froot
froot=( `read_file_column $input_table 2 1` )
pass=`echo ${froot[@]}`
cnt=( `array_size "$pass"` )

for ((i=0;i<$cnt;++i)); do
    if [ $(find_string ${froot[${i}]} ${2}) -eq 1 ]; then
	break
    fi
done

# Run the DAP on this file
if [ -e ${3} ]; then
    echo "WARNING: Overwriting existing spool file: ${3}"
    rm -f ${3}
fi
cmd='idl -quiet -e "manga_dap, $i, /nolog" > $3'
eval $cmd

# End time
date

# END
####################################################################################################
exit 0


