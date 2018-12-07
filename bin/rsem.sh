#!/bin/bash -x
set -eu -o pipefail
bam=$1 #star_align/$SID/Aligned.toTranscriptome.out.bam"
gdir=$2
threads=$3
paired=$4 #0 is single and 1 paired

SID=$(basename $(dirname $bam))
pairopt=""
if (( $paired == 1 )); then
    pairopt="--paired-end"
fi
#we have added the version that can work with single ends data as well
rsem-calculate-expression \
    $pairopt \
    -p $threads\
    --no-bam-output\
    --forward-prob 0.5\
    --seed 12345\
    --bam $bam\
    $gdir/rsem_index/genome\
    rsem/$SID
