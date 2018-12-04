#!/bin/bash
set -euo pipefail
threads=$1
bowtie2-build --threads $threads genome.fa bowtie2_index/genome
