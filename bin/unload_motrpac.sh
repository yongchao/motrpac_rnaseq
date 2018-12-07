#!/bin/bash 
set -eu -o pipefail
#remove the environmetnal variables from the PATH
#Running this command itself is not going to change the  environmetnal variables
#you have to run
#export $(unlod_motrpac.sh)
remove_path(){
    p=$2
    echo $1 | awk -v p=$p 'BEGIN{RS=ORS=":"}; $0!=p'|
	awk '{sub(":$","",$0);print}'
}
set +u 
if [ "${MOTRPAC_root}x" != x ]; then
    PATH=$(remove_path $PATH $MOTRPAC_root/bin)
    PATH=$(remove_path $PATH $MOTRPAC_conda/python3/bin)
    PATH=$(remove_path $PATH $MOTRPAC_conda/python2/bin)
fi
set -u

echo export PATH=$PATH
echo unset MOTRPAC_root
echo unset MOTRPAC_conda
echo unset MOTRPAC_refdata

