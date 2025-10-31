#!/bin/bash 
set -eux -o pipefail

set +x
ml bcl2fastq/2.20.0.422 #the  plain old versuion doesn't work
set -x

mkdir -p interop ##for storing the file IndexMetricsOut.bin,
##need to check if it's working

mkdir -p log

threads=6
mem=$((32000/threads+4000))

bsub -K -q premium -n $threads -R "rusage[mem=$mem] span[hosts=1]" \
     -W 24:00 -P acc_sealfs01a\
      -oo log/bcl2fastq.out  -eo log/bcl2fastq.err -J bcl2fastq \
    bcl2fastq --sample-sheet ../SampleSheet.csv \
    -p $threads \
    --use-bases-mask Y*,I8Y*,I8N*,Y* \
    -R  ../bcl_rawdata \
    --mask-short-adapter-reads 0 \
    --minimum-trimmed-read-length 0 \
    -o . \
    --interop-dir interop "$@"

ln -fs ./Reports/html/*/all/all/all/laneBarcode.html .
