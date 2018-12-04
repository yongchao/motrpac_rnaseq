#!/bin/bash -x
set -eu -o pipefail

gdir=$1
threads=$2
shift 2
#The remaining one or two parameters are the input fastq files (can be one or two)
SID=$(basename $1 _R1.fastq.gz)

mkdir -p star_align/$SID

#why is this required
ulimit -v 41000000000

STAR  --genomeDir $gdir/star_index\
      --sjdbOverhang  99\
      --readFilesIn "$@"\
      --outFileNamePrefix star_align/${SID}/\
      --readFilesCommand zcat \
      --outSAMattributes NH HI AS NM MD nM\
      --outFilterType BySJout\
      --runThreadN $threads\
      --outSAMtype BAM SortedByCoordinate\
      --quantMode TranscriptomeSAM
      #--genomeFastaFiles $gdir/genome.fa\
      #--sjdbGTFfile $gdir/genome.gtf\
      #--twopassMode Basic
      #The above three options not working for the memory

      #--twopassMode Basic
      #--limitBAMsortRAM 15000000000
      #--genomeLoad NoSharedMemory\ doesn't seem to work, but it's the default anyway
      #--genomeLoad LoadAndKeep \ the MOP is wrong here
