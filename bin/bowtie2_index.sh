#!/bin/bash
set -euo pipefail
threads=$1
log=$2
bowtie2-build genome.fa bowtie2_index/genome >&{log}
