#!/bin/bash -x
set -eu -o pipefail 
SID=$1
gdir=$2
threads=$3

if [[ -e fastq/${SID}_R2.fastq.gz ]]; then
    cmd="-1 fastq/${SID}_R1.fastq.gz -2 fastq/${SID}_R2.fastq.gz"
else
    cmd="-U fastq/${SID}_R1.fastq.gz"
fi

#module load bowtie2/2.2.8, using the conda version 2.2.4
set +x
module load samtools/1.4.1
set -x

genome=$(basename $gdir)
genome=${genome%%_*}
##remove the version

genome=${genome%%[0-9]*}
rRNA=/sc/orga/projects/sealfs01a/motrpac/refdata/rRNA/${genome}_rRNA

bowtie2 --local -p $threads -s --sam-nohead -x $rRNA $cmd >rRNA/$SID.sam 2>rRNA/$SID.txt

##The following is only needed when when want to investigate the bam file for details
cd rRNA
samtools view -b $SID.sam -o ${SID}_unsorted.bam

tmpdir=${SID}_samtools_tmp
mkdir -p $tmpdir

samtools sort --threads $threads -T $tmpdir -o $SID.bam ${SID}_unsorted.bam
samtools index ${SID}.bam

##clean out the intermediate files
rm $SID.sam
rm ${SID}_unsorted.bam
rm -rf samtools_sortq

