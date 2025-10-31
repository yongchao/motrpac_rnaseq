#!/bin/bash -x
set -eu -o pipefail
bam=$1 #star_align/{sample}/Aligned.sortedByCoord.out.bam
gdir=$2
threads=$3
paired=$4 #0 is single and 1 paired
tmpdir_root=$5 #the temporary file folder
SID=$(basename $(dirname $bam))
pairopt=""
if (( $paired == 1 )); then
    pairopt="-p"
fi

tmpdir=$(mktemp -d -p $tmpdir_root featureCounts.${SID}.XXX)  

featureCounts -T $threads --tmpDir $tmpdir  -a $gdir/genome.gtf -o featureCounts/$SID $pairopt -M --fraction $bam
rm -rf $tmpdir
