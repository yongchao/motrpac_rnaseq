#!/bin/bash
set -eu -o pipefail
SID=$1
grp=$2
gdir=$3
module load java/1.8.0_111
module load picard/2.7.1
strand=SECOND_READ_TRANSCRIPTION_STRAND #The MOP specify FIRST_READ_TRANSCRIPTION_STRAND, this depends on libary prep
java -jar $PICARD CollectRnaSeqMetrics\
     I=star_align/${SID}/Aligned.sortedByCoord.out.bam \
     O=qc53/${SID}${grp}.RNA_Metrics\
     MINIMUM_LENGTH=50 \
     REF_FLAT=$gdir/qc53_ref/genome_flat${grp}.txt\
     STRAND_SPECIFICITY=$strand \
     RRNA_FRAGMENT_PERCENTAGE=0.3 \
     >& qc53/log/${SID}${grp}.log
