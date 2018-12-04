#!/bin/bash -x
set -euo pipefail
fqfile=$1 #possily from fastq (or fastq_trim), fastq_raw
odir=$2 #possiblly fastqc, or fastqc_raw
fastqc -o $odir $fqfile 
