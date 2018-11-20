#!/bin/bash 
set -eu -o pipefail

#remove a path from the PATH
remove_path(){
    p=$2
    echo $1 | awk -v p=$p 'BEGIN{RS=ORS=":"}; $0!=p'|
	awk '{sub(":$","",$0);print}'
}
set +u 
if [ "${MOTRPAC_ROOT}x" != x ]; then
    PATH=$(remove_path $PATH $MOTRPAC_ROOT/bin)
fi
if [ "${CONDA_ROOT}" != x ]; then
    PATH=$(remove_path $PATH $CONDA_ROOT/bin)
fi
set -u

PATH=$(remove_path $PATH /sc/orga/projects/sealfs01a/motrpac/bin)
PATH=$(remove_path $PATH /sc/orga/projects/sealfs01a/conda/python3/bin)

echo export PATH=$PATH
echo unset MOTRPAC_ROOT
echo unset CONDA_ROOT
echo unset MOTRPAC_TMP


