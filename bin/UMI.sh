#!/bin/bash -x
set -euo pipefail
#output the de-duplicated for NuGEN UMI stuffs
#Default is 8 bases, this can be changed as well
bam=$(readlink -e $1)
trans_bam=$(readlink -e $2)
paired=$3
pairopt=""
if (( $paired == 1 )); then
    pairopt="-2"
fi

SID=$(basename $(dirname $bam))
set +e #head has a problem with this
len=$(samtools view $bam |head -1 |awk '{umi=gensub("^.*:","","a",$1); print length(umi)}')
set -e

#filter uniquelly mapped reads

#The default was unitelligently set to 6
cd star_align/UMI_$SID

tmpdir=$(readlink -e $MOTRPAC_ROOT/tmpdir) #this may avoid the nudup.py named pipe problems
bam_uni=${SID}_uniq.bam

bam_uniq.sh $bam $bam_uni $paired 
python2 $MOTRPAC_ROOT/nugen/nudup.py $pairopt -s $len -l $len --rmdup-only -o $SID -T $tmpdir $bam_uni
mv $SID.sorted.dedup.bam Aligned.sortedByCoord.out.bam

#This file needs to be sorted first
samtools sort -o trans_sorted.bam $trans_bam
bam_uniq.sh trans_sorted.bam $bam_uni $paired
python2 $MOTRPAC_ROOT/nugen/nudup.py $pairopt -s $len -l $len --rmdup-only -o ${SID}_trans -T $tmpdir $bam_uni

#clean up
rm trans_sorted.bam $bam_uni
#rsem requires names to be sorted
samtools sort -n -o Aligned.toTranscriptome.out.bam  ${SID}_trans.sorted.dedup.bam 
rm ${SID}_trans.sorted.dedup.bam 


