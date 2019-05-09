#!/bin/bash -x
set -eu -o pipefail
bam=$1 #star_align/$SID/Aligned.toTranscriptome.out.bam"
gdir=$2
threads=$3
paired=$4 #0 is single and 1 paired
tmpdir_root=$5 #tmp folder
SID=$(basename $(dirname $bam))
pairopt=""
if (( $paired == 1 )); then
    pairopt="--paired-end"
fi
#we have added the version that can work with single ends data as well
tmpdir=$(mktemp -d -p $tmpdir_root rsem.${SID}.XXX)

rsem-calculate-expression \
    $pairopt \
    -p $threads\
    --no-bam-output\
    --forward-prob 0.5\
    --seed 12345\
    --bam $bam\
    --temporary-folder $tmpdir\
    $gdir/rsem_index/genome\
    rsem/$SID
rm -rf $tmpdir
