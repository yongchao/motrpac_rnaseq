#!/bin/bash
set -euo pipefail
bam=$1
bam_u=$2
paired=$3 #this could infer from the data, but tedious
if (( paired ==1)); then
    samtools view -F 0xF0c -f 0x03 $bam -o $bam_u
else
    samtools view -F 0xF05 $bam -o $bam_u
fi
#uniquely mapped flags: to do list for both rsem and star_align
#-F 0x100+0x200+0x400+0x800+0x004+0x008=0xF0c and -f 0x03, only works for paired
#-F 0x100+0x200+0x400+0x800+0x004+0x001=0xF05 #for single
