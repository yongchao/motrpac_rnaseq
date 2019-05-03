#!/bin/bash
set -euo pipefail

#Note the input only includes R1 or R2 fastq files from the fastq_raw

#index_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
#univ_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#There adapters are too long, 3'end usually has lower quality

#Following trim_galore approach, taking the first 13 bases
#Maybe we should go with trim_galore

index_adapter=AGATCGGAAGAGC
univ_adapter=AGATCGGAAGAGC

#The following scrips make it possible to supplie the adpaters for non-mtropac projects
while getopts a:A:h o 
do      
    case "$o" in
        a) index_adapter="$OPTARG";;
	A) univ_adapter="$OPTARG";;
        h) echo "Usage: $0 [-h] [-a index_adapter ] [-A univ_adapter ] fq1 [fq2]"
           echo '-h: print help'
	   echo '-a index_adapater: set the index_adapter (default=AGATCGGAAGAGC)'
	   echo '-A universal_adapater: set the univ_adapter (default=AGATCGGAAGAGC)'
           exit 0;;
    esac
done

shift $((OPTIND-1))

#$2, when present, must be the same as $1 with R1 replaced by R2, not enforced here for now
R1=$(basename $1)
R2=""
index2=""
cmdo2=""
cmd_tooshort2=""

if (( $# == 2 )); then
    R2=$(basename $2)
    index2="-A $univ_adapter"
    cmdo2="-p fastq_trim/$R2"
    cmd_tooshort2="--too-short-paired-output fastq_trim/tooshort/$R2"
fi
mkdir -p fastq_trim/tooshort    
cutadapt \
    -a $index_adapter $index2 \
    -o fastq_trim/$R1 $cmdo2\
    -m 20 \
    --too-short-output fastq_trim/tooshort/$R1 $cmd_tooshort2\
    "$@"


    
   
