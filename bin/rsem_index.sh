#!/bin/bash
set -euo pipefail
log=$1
#As of rsem 1.3.1. There is no support for the multi threads. The threads
#are only supported for the optionally star and bowtie index, which is
#not needed here
rsem-prepare-reference --gtf genome.gtf genome.fa rsem_index/genome >& &log
