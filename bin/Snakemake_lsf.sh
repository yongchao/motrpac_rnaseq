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
while getopts j:W:P:q:n:R:m:h o 
do      
    case "$o" in
        j) jsonfile="$OPTARG";;
        W) optW="$OPTARG";;
	P) optP="$OPTARG";;
	q) optq="$OPTARG";;
	n) optn="$OPTARG";;
	R) optR="$OPTARG";;
	m) cmdm="-m $OPTARG";;
        h) echo "Usage: $0 [-h] [-j json] [-W wall_time] [-P account] [-q queue [-p proj] [-- -s snakefile and other snakemake options]"
           echo '-h: print help'
	   echo "-j: the json file (default $MOTRPAC_ROOT/config/lsf.json)"
           echo '-W: the wall time (default 4:00)'
	   echo '-P: the account (default acc_sealfs01a)'
	   echo '-q: the queue (default premimum)'
	   echo '-n: the number of CPUs (default 1)'
	   echo '-R: resources (default \"rusage[mem=4000]\")'
	   echo '-m: restrict the machine types (mothra, mandra,bode)'
           exit 0;;
    esac
done

shift $((OPTIND-1))

mkdir -p log/cluster

set -x

snakemake -j 999 --cluster-config $jsonfile --cluster "bsub $cmdm -W $optW -P $optP -q $optq -n $optn -R $optR -oo {cluster.output} -eo {cluster.error} -J {cluster.name}" "$@"
