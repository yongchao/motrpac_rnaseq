#!/bin/bash 
set -eu -o pipefail

#remove the duplicate term with only one variable (that do not contain space or other special characters)
de_duplicate(){
    echo $1| awk 'BEGIN{RS=ORS=":"};!a[$1]++'|
	awk '{sub(":$","",$0);print}'
}
#export -f de_duplicate

optr=/sc/orga/projects/sealfs01a/motrpac
optc=/sc/orga/projects/sealfs01a/conda/python3

while getopts r:h o 
do      
    case "$o" in
	h) echo "Usage: $0 [-h] [-r motrpac_code_root] [-c conda]"
	   echo '-h: print help'
	   echo '-r: specify the root of the code base for motropac (default /sc/orga/projects/sealfs01a/motrpac for production)'
	   echo '    it can be your home folder eg ~/github/sealfonlab/motrpac for testing'
	   echo '-c: specify the conda folder (defualt /sc/orga/projects/sealfs01a/conda/python3 for production)'
	   exit 0;;
	r) optr="$OPTARG";;
	c) optc="$OPTARG";;
    esac
done

##Simple testing the $optr is pointing the right motrpac root folder
if [ ! -d $optr/refdata ]; then
    echo "$optr can not be used to specify as the root for the motrpac code as it didn't have refdata subfolder"
    exit 1
fi

if [ ! -d $optc/conda-meta ]; then
    echo "$optc can not be used to specify as the conda environment as it didn't have conda-meta subfolder"
    exit 1
fi

echo export MOTRPAC_ROOT=$optr
echo export CONDA_ROOT=$optc
echo export PATH=$(de_duplicate $optc/bin:$optr/bin:$PATH)
