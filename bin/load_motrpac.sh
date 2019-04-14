#!/bin/bash 
set -eu -o pipefail
# Running this script WILL NOT change the enviromental variables

# You need to run the comand below to export the the environmental variables $MOTRPAC_root, $MOTRPAC_conda,$MOTRPAC_refdata,$PATH
# first try the command bin/load_motrpac.sh [-c conda_dir] [-r refdata_dir] to identify problems, 
# and then running the command below to export the variables, or you can just copy and paste the output from the above command to the terminal
# export $(bin/load_motrpac.sh [-c conda_dir] [-r refdata_dir])

#load_motrpac.sh might be a link
root=$(dirname $(readlink -m $(dirname $0))) #find out the real motrpac root folder
conda=/sc/hydra/projects/sealfs01a/conda #sinai default
refdata=/sc/hydra/projects/sealfs01a/motrpac_refdata #sinai default
while getopts hc:r: o 
do      
    case "$o" in
	c) conda="$OPTARG";;
	r) refdata="$OPTARG";;
	h) echo "Usage: $0 [-h] [-c conda_dir] [-r refdata_dir] "
	   echo '  : prints out the environment variables that needs to be exported and '
	   echo '    the folder structure is correct'     
	   echo '-h: print help'
	   echo '-c conda_dir: setting the folder for the conda (default for sinai setup)'
	   echo '-f refdata_dir: setting the folder for the genome references (default for sinai setup)'
	   exit 0;;
    esac
done
if [[ ! -d $refdata  ]]; then
    echo "The refdata folder is not a folder"
    exit 1
fi

if [[ ( ! -d $root/config || ! -d $root/nugen ) || (! -f $root/bin/unload_motrpac.sh || ! -f $root/sample_sub.snakefile) ]]; then
    echo "The motrpac root folder does not have the right folder structure"
    exit 1
fi

if [[ ! -d $conda/python2/conda-meta || ! -d $conda/python3/conda-meta ]]; then
    echo "The conda folder does not have the right folder structure"
    exit 1
fi
#remove the duplicate paths (assuming no path contain space or other special characters)
de_duplicate(){
    echo $1| awk 'BEGIN{RS=ORS=":"};!a[$1]++'|
	awk '{sub(":$","",$0);print}'
}

#We need to unload the previous enviroment setup
#This segment is almost the same as script unload_motrpac.sh
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
##########

echo export PATH=$(de_duplicate $root/bin:$conda/python3/bin:$conda/python2/bin:$PATH)
echo export MOTRPAC_root=$root
echo export MOTRPAC_conda=$conda
echo export MOTRPAC_refdata=$refdata
