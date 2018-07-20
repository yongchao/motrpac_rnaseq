#!/bin/bash -x 
prep_sample.py  | sort_skip.sh 1 -k3,3n -k1,1 > ../sample_info.txt
cd ..

g=hg38_gencode_v28
if [[ "$#" -ge  1 ]]; then 
    g=$1
fi

printf "genomedir:\n    $MOTRPAC_ROOT/refdata/$g\n" >config.yaml

mkdir -p log/cluster

