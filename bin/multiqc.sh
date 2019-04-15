#!/bin/bash
set -eu -o pipefail
name=$1
shift 1
#The output is always in the multiqc folder
#The command is used twice, once in rule pre_align_QC, and another time in rule post_align_QC

#the following setup is necessary
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export TMPDIR=tmpdir
mkdir -p $TMPDIR

multiqc -d \
	-f \
	-n $name \
	-o multiqc \
	"$@"

