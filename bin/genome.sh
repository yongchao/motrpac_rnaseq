#!/bin/bash
set -euo pipefail
#we are currently at the root of the genome folder
#generate the genome.fa and genome.gtf from the source folder
#This is brute force, though it doesn't not change for gencode data
zcat source/*.fa.gz | fixchr4ensembl.sh >genome.fa
zcat source/*.gtf.gz | fixchr4ensembl.sh -g |sort -k1,1 -k4,4n -k5,5n>genome.gtf
