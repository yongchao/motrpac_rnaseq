#!/bin/bash -x
set -eu -o pipefail
bam=$1 #star_align/{sample}/Aligned.sortedByCoord.out.bam
gdir=$2
threads=$3
paired=$4 #0 is single and 1 paired
tmpdir_root=$5 #the temporary file folder
SID=$(basename $(dirname $bam))
pairopt=""
if (( $paired == 1 )); then
    pairopt="-p"
fi

tmpdir=$(mktemp -d -p $tmpdir_root featureCounts.${SID}.XXX)

sam=$tmpdir/$SID.sam
sam_sorted=$tmpdir/${SID}_sorted.sam
sam_full=$tmpdir/${SID}_full.sam
		     
#samtools sort -m 5.5G -n -T $tmpdir/samsort -@ $threads -o $bam_paired $bam
#it generates too many small files and it is not working well

samtools view -@ $threads $bam >$sam

export LC_ALL=C; export LC_LANG=C
sort $sam -k1,1 -o $sam_sorted -T $tmpdir --parallel=$threads -S 5G
samtools view -H $bam >$sam_full
cat $sam_sorted >>$sam_full

#this seems promising

featureCounts -T $threads --tmpDir $tmpdir  -a $gdir/genome.gtf -o featureCounts/$SID $pairopt -M --fraction $sam_full --donotsort
rm -rf $tmpdir
