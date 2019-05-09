#!/bin/bash -x
set -eu -o pipefail

gdir=$1
threads=$2
tmpdir_root=$3
shift 3
#The remaining one or two parameters are the input fastq files (can be one or two)
SID=$(basename $1 _R1.fastq.gz)

mkdir -p star_align/$SID

tmpdir=$(mktemp -d -p $tmpdir_root star_align.${SID}.XXX)/tmp 
#why is this required
#ulimit -v 41000000000

STAR  --genomeDir $gdir/star_index\
      --sjdbOverhang  100\
      --readFilesIn "$@"\
      --outFileNamePrefix star_align/${SID}/\
      --readFilesCommand zcat \
      --outSAMattributes NH HI AS NM MD nM\
      --runThreadN $threads\
      --outSAMtype BAM SortedByCoordinate\
      --outFilterType BySJout\
      --quantMode TranscriptomeSAM\
      --outTmpDir $tmpdir
rm -rf $tmpdir
#--outFilterType BySJout\ with the the problem downstream analysis for the UMI as the paired reads have been removed
