#!/bin/bash 
set -euo pipefail

#This script is better to run interactively to find where the problem is.

#Please note that SampleSheet.csv needs to have the adpater setting removed
#as we are going to use cutadapt to remove the adapter anyway.

mkdir bcl2fastq
cd bcl2fastq

# $bcl_rawdata is defined to a folder that contains the raw data of the bcl run folder
# It can be defined here by
# export bcl_raw=/path/to/bcl_raw_data

bcl2fastq --sample-sheet SampleSheet.csv \
	  --use-bases-mask Y*,I8Y*,I*,Y* \
	  -R  $bcl_rawdata \
	  -o . --minimum-trimmed-read-length 0 \
	  --mask-short-adapter-reads 0 

#The following is necessary to rename R2 to I1 and R3 to R2 files
find . -name "*_R2_001.fastq.gz" |awk '{print "mv "$0" "gensub("_R2_001.fastq.gz$","_I1_001.fastq.gz",1)}' |bash -x
find . -name "*_R3_001.fastq.gz" |awk '{print "mv "$0" "gensub("_R3_001.fastq.gz$","_R2_001.fastq.gz",1)}' |bash -x

#If we only have a single lane, we can build soft links as below, otherwise we need to merge the fastq files from different lanes.

#In Illumina, the fastq file has the format:
#SID_S[0-9]+_L00[1-8]_[IR][12]_001.fastq.gz, we assume SID is unique for different samples, which can be done by changing file SampleSheet.csv
cd ..
mkdir fastq_raw
cd fastq_raw
find ../bcl2fastq -name "*_001.fastq.gz" |awk 'BEGIN{FS="/"};{print "ln -s "$0" "gensub(/_S[0-9]+_L00[1-8](_[IR][12])_001/,"\\1",1,$(NF))}' |bash -x
