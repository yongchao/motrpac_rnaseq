#!/bin/bash 
set -eu -o pipefail

#Runing this can not change the enviromental variables
#The was to export the the enviromental variables at the MOTRPAC_ROOT folder
#export $(bin/load_motrpac.sh -c $conda/python3)

optr=$(dirname($(dirname $(readlink -m $0)))) #find out the real motrpac root folder
optc=/sc/orga/projects/sealfs01a/conda/python3 #This is Sinai's default
optt=/sc/orga/scratch/$USER/tmp #This is Sinai's tmp folder

while getopts r:ht: o 
do      
    case "$o" in
	h) echo "Usage: $0 [-h] [-r motrpac_code_root] [-c conda]"
	   echo '-h: print help'
	   echo '-r: specify the root of the code base for motropac (default is the github code root of motrpac)'
	   echo '    This can be changed to some frozen code base'
	   echo '-c: specify the conda folder (default /sc/orga/projects/sealfs01a/conda/python3)'
	   echo '-t: specify the temp folder (default /sc/orga/scratch/$USER/tmp)'
	   exit 0;;
	r) optr="$OPTARG";;
	c) optc="$OPTARG";;
	t) optt="$OPTARG";;
    esac
done

##Simple testing the $optr to make sure it points to the right motrpac root folder
if [ ! -d $optr/refdata ]; then
    echo "$optr can not be used to specify as the root for the motrpac code as it didn't have refdata subfolder"
    exit 1
fi

if [ ! -d $optc/conda-meta ]; then
    echo "$optc can not be used to specify as the conda environment as it didn't have conda-meta subfolder"
    exit 1
fi

mkdir -p $optt

echo export MOTRPAC_ROOT=$optr
echo export CONDA_ROOT=$optc
echo export MOTRPAC_TMP=$optt

#remove the duplicate paths (assuming no path contain space or other special characters)
de_duplicate(){
    echo $1| awk 'BEGIN{RS=ORS=":"};!a[$1]++'|
	awk '{sub(":$","",$0);print}'
}

echo export PATH=$(de_duplicate $optc/bin:$optr/bin:$PATH)
