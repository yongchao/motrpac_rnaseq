#!/bin/bash
set -eu -o pipefail
#if -b set, then computing the bam file in the current folder, otherwise just the qc metrics
#if -d dir set, then all of the results including bam file is saved at the folder dir
optb=0
dir="."
while getopts bhd: o 
do      
    case "$o" in
        b) optb=1;;
	d) dir="$OPTARG";;
        h) echo "Usage: $0 [-h] [-b [-d dir] ] gref threads fq1 [fq2] "
           echo '-h: print help'
	   echo '-b: compute the bam files'
	   echo '-d dir: the bam file is saved dir'
           exit 0;;
    esac
done

shift $((OPTIND-1))

gref=$1 #bowite2 ref
threads=$2
shift 2

#The remaining one or two parameters are the input fastq files (can be one or two)
SID=$(basename $1 _R1.fastq.gz)
if (($# == 2)); then #paired
    cmd="-1 $1 -2 $2"
else
    cmd="-U $1"
fi

sam=$dir/$SID.sam
#Using local in case the the adapter is not working
bowtie2 -p $threads $cmd -x $gref --local -S $sam

sleep 5 #sometimes it takes time for the cluster to update
#sometimes $sam is not generated due to zero alignment
if [[ -f $sam && -s $sam ]]; then
    if  (($optb==1)); then
	bam_unsorted=$dir/${SID}_unsorted.bam
	bam=$dir/$SID.bam
	samtools view -@ {threads} -b $sam -o $bam_unsorted
	samtools sort -@ {threads} -o $bam $bam_unsorted
	samtools index $bam
	rm $bam_unsorted
    fi
    rm $sam
fi
