#!/bin/bash
set -euo pipefail
##Assuming genome.fa and genome.gtf is in the current folder
threads=$1
STAR --runThreadN $threads\
     --runMode genomeGenerate\
     --genomeDir star_index\
     --genomeFastaFiles genome.fa\
     --sjdbGTFfile genome.gtf\
     --sjdbOverhang 100\
     --outFileNamePrefix star_index/ 
