#!/bin/bash
set -euo pipefail
##Assuming genome.fa and genome.gtf is in the current folder
threads=$1
log=$2
STAR --runThreadN $threads\
     --runMode genomeGenerate\
     --genomeDir star_index\
     --genomeFastaFiles genome.fa\
     --sjdbGTFfile genome.gtf\
     --sjdbOverhang 99\
     --outFileNamePrefix star_index/ >& $log
