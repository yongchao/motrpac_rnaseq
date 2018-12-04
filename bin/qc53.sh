#!/bin/bash
set -eu -o pipefail
bam=$1 #star_align/{sample}_Aligned.sortedByCoord.out.bam
gdir=$2
SID=$(basename $(dirname $bam))
strand=FIRST_READ_TRANSCRIPTION_STRAND #The MOP specify FIRST_READ_TRANSCRIPTION_STRAND
picard CollectRnaSeqMetrics\
     I=$bam \
     O=qc53/${SID}.RNA_Metrics\
     MINIMUM_LENGTH=50 \
     REF_FLAT=$gdir/qc53_ref/genome_flat.txt\
     STRAND_SPECIFICITY=$strand \
     RRNA_FRAGMENT_PERCENTAGE=0.3 \
     >& qc53/log/${SID}.log
