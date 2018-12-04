#!/bin/bash 
set -eu -o pipefail

#Before running this, makesure refdata, conda and tmpdir under the MOTRPAC_ROOT have already been soft-linked to the right folders

#Running this WILL NOT change the enviromental variables

#You need to run the comand below to export the the enviromental variables at the $MOTRPAC_ROOT and $PATH
#export $(bin/load_motrpac.sh)

root=$(dirname $(dirname $(readlink -m $0))) #find out the real motrpac root folder

while getopts h o 
do      
    case "$o" in
	h) echo "Usage: $0 [-h] "
	   echo '-h: print help'
	   echo '  : prints out the environment variables that needs to be exported and '
	   echo '    the folder structure is correct'
	   exit 0;;
    esac
done

##Simple testing the $root to make sure it points to the right motrpac root folder
#to makre sure refdata is not pointing a wrong folder -d $root/refdata/globin
if [ -L $root/redata ]; then
    echo "MOTRPAC_ROOT requires a softlink to refdata"
    exit 1
fi

if [[ ! -L $root/conda || ! -d $root/conda/python2/conda-meta || ! -d $root/conda/python3/conda-meta ]]; then
    echo "MOTRPAC_ROOT requires a softlink to conda with the right folder structure"
    exit 1
fi

if [[ ! -L $root/tmpdir ]]; then
    echo "MOTRPAC_ROOT requires softlinke tmpdir"
    exit 1
else
    #find out the space storage
    space=$(df --direct $root/tmpdir |tail -1|awk '{print $4}')
    if (( space < 1000000000 )); then
	echo "MOTRPAC_ROOT/tmpdir needs to be at least 100G space"
	exit 1
    fi
fi

#remove the duplicate paths (assuming no path contain space or other special characters)
de_duplicate(){
    echo $1| awk 'BEGIN{RS=ORS=":"};!a[$1]++'|
	awk '{sub(":$","",$0);print}'
}

echo export PATH=$(de_duplicate $root/bin:$root/conda/python3/bin:$root/conda/python2/bin:$PATH)
echo export MOTRPAC_ROOT=$root
