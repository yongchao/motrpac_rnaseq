#!/bin/bash
set -euo pipefail
log=$1
module load ucsc-utils/2015-04-07
cd qc53_ref
gtfToGenePred -ignoreGroupsWithoutExons -genePredExt ../genome.gtf genome_pred.txt
awk 'BEGIN{FS=OFS="\t"};{print $12, $1, $2,$3,$4,$5,$6,$7,$8,$9,$10}' genome_pred.txt>genome_flat.txt
rm genome_pred.txt
