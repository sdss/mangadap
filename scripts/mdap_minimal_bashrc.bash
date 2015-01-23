export MDAP_DEV=/Users/andrewsb/manga/mangadap/brett_v0.93

export IDLUTILS_DIR=/Users/andrewsb/.idl/idlutils
export IDL_DIR=/Applications/itt/idl/idl81
export IDL_PATH=$IDL_DIR/lib
export PATH=$IDLUTILS_DIR/bin:$PATH

if [ -e $HOME/.idlstartup ]; then
  export IDL_STARTUP=$HOME/.idlstartup
fi