#!/bin/bash
set -eu -o pipefail
bam=$1 #star_align/{sample}_Aligned.sortedByCoord.out.bam
SID=$(basename $(dirname $bam))
PICARD=$(dirname $(readlink -e $(which picard)))/picard.jar
#In future try to obtain from the resources
java -Xmx35000M -jar $PICARD MarkDuplicates \
     I=$bam \
     O=mark_dup/${SID}_markedDup.bam \
     CREATE_INDEX=true \
     VALIDATION_STRINGENCY=SILENT \
     ASSUME_SORT_ORDER=coordinate\
     M=mark_dup/${SID}.dup_metrics \
     REMOVE_DUPLICATES=false 
rm mark_dup/${SID}_markedDup.bam*
