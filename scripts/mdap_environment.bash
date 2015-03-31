#!/bin/bash

# This is the path to the equivalent of 
#   $MANGAWORK_DIR/manga/spectro
# at Utah, with your local copy of the DRP products
export LOCAL_MANGA_PATH=/Users/westfall/Work/MaNGA

export MANGA_SPECTRO_REDUX=$LOCAL_MANGA_PATH/redux
export MANGA_SPECTRO_ANALYSIS=$LOCAL_MANGA_PATH/analysis
export MANGA_SPECTRO_ANALYSIS_PLOTS=$LOCAL_MANGA_PATH/analysis/plots

#export MANGADRP_VER=v1_0_0
#export MANGADRP_VER=v1_1_2
export MANGADRP_VER=v1_2_0
#export MANGADAP_VER=v0_9_0
export MANGADAP_VER=trunk_v1_2_0_mpl3beta

export MANGADAP_DIR=$MDAP_DEV

return 0

