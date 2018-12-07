#!/bin/bash
set -euo pipefail
bam=$1  #star_align/$SID/Aligned.sortedByCoord.out.bam
SID=$(basename $(dirname $bam))
bam2=star_align/$SID/${SID}_primary.bam
#Assuming bam is sorted
#filter out the primary alignment
samtools view -b -F 0x900 $bam -o $bam2
samtools index $bam2
samtools idxstats $bam2 > star_align/$SID/chr_info.txt
rm $bam2 $bam2.bai
