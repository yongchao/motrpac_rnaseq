#!/bin/bash 
set -euo pipefail
IFS=$'\n\t'
if (( $# >0 )); then
    n=$1
else
    n=0
fi

if (( $# >1 )); then
    k=$2	
else
    k=0
fi
cut -f $n | tail -n +$((1+$k))
