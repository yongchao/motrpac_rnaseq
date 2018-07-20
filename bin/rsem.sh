#!/bin/bash -x
set -eu -o pipefail
SID=$1
gdir=$2
threads=$3
 #module load rsem/1.2.29, using the conda version
pairopt=""
if [ -e fastq/${SID}_R2.fastq.gz ]; then
    pairopt="--paired-end"
fi
#we have added the version that can work with single ends data as well
rsem-calculate-expression \
    $pairopt \
    -p $threads\
    --no-bam-output\
    --forward-prob 0.5\
    --seed 12345\
    --bam  star_align/$SID/Aligned.toTranscriptome.out.bam\
    $gdir/rsem_index/genome\
    rsem/$SID
