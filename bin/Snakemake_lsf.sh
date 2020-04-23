#!/bin/bash 
set -eux -o pipefail

#-j jsonfile
#-W wall time
#-P project
#-q alloc
#-n nCPU
#-R resource
#-m mothra
jsonfile=$MOTRPAC_root/config/lsf.json
optn="{cluster.nCPUs}"
optW="{cluster.time}"
optP="{cluster.account}"
optq="{cluster.queue}"
optR="{cluster.resources}"
cmdm=""
cmdx=""
cmdK=""

while getopts j:W:P:q:n:R:m:hxK o 
do      
    case "$o" in
        j) jsonfile="$OPTARG";;
        W) optW="$OPTARG";;
	P) optP="$OPTARG";;
	q) optq="$OPTARG";;
	n) optn="$OPTARG";;
	R) optR="$OPTARG";;
	m) cmdm="-m $OPTARG";;
	x) cmdx="-sla Sealfon_sla";;
	K) cmdK="-K";;
        h) echo "Usage: $0 [-h] [-j json] [-W wall_time] [-P account] [-K] [-q queue] [-p proj] [-- -s snakefile and other snakemake options]"
           echo '-h: print help'
	   echo "-j: the json file (default $MOTRPAC_ROOT/config/lsf.json)"
           echo '-W: the wall time (default {cluster.time})'
	   echo '-P: the account (default acc_sealfs01a)'
	   echo '-q: the queue (default premimum)'
	   echo '-n: the number of CPUs (default {cluster.nCPUs})'
	   echo '-R: resources (default {cluster.resources})'
	   echo '-m: restrict the machine types (mothra, mandra,bode)'
	   echo '-K: no exit'
           exit 0;;
    esac
done

shift $((OPTIND-1))

mkdir -p log/cluster

set -x

snakemake -j 400 --cluster-config $jsonfile --cluster "bsub $cmdm  -W $optW -P $optP $cmdx $cmdK -q $optq -n $optn -R $optR -oo {cluster.output} -eo {cluster.error} -J {cluster.name}" "$@"

