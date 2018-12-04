#!/bin/bash
set -euo pipefail

#Note the input only includes R1 or R2 fastq files from the fastq_raw

#index_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# univ_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#There adapters are too long, 3'end usually has lower quality

#Following trim_galore approach, taking the first 13 bases
#Maybe we should go with trim_galore

index_adapter=AGATCGGAAGAGC
 univ_adapter=AGATCGGAAGAGC
#Future plan is to make the adpaters can be supplied on the fly
#$2, when present, must be the same as $1 with R1 replaced by R2, not enforced here for now
R1=$(basename $1)
R2=""
cmdo2=""
cmd_tooshort2=""

if (( $# == 2 )); then
    R2=$(basename $2)
    cmdo2="-p fastq_trim/$R2"
    cmd_tooshort2="--too-short-paired-output fastq_trim/tooshort/$R2"
fi
mkdir -p fastq_trim/tooshort    
cutadapt \
    -a $index_adapter \
    -A $univ_adapter \
    -o fastq_trim/$R1 $cmdo2\
    -m 20 \
    --too-short-output fastq_trim/tooshort/$R1 $cmd_tooshort2\
    "$@"


    
   
