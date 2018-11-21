#!/bin/bash
set -euo pipefail
#we are currently at the root of the genome folder
#generate the genome.fa and genome.gtf from the source folder
#This is brute force, though it doesn't not change for gencode data
zcat *.fa.gz | source/fixchr4ensembl.sh >genome.fa
zcat *.gtf.gz | source/fixchr4ensembl.sh |sort -k1,1 -k2,2n -k3,3n genome.gtf
