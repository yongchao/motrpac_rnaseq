#!/bin/bash
#-g for gtf file
#otherwise for fa file
optG=0
while getopts g o
do
    case "$o" in
	g) optG=1;;
    esac
done

shift $((OPTIND-1))
if [[ $optG -eq 1 ]]; then
    sed 's/^\([0-9]\+\)/chr\1/' |
    sed 's/^MT/chrM/'|
    sed 's/^\([XY]\)/chr\1/' 
else
#for the fasta file
    sed 's/^>\([0-9]\+\)/>chr\1/' |
    sed 's/^>MT/>chrM/'|
    sed 's/^>\([XY]\)/>chr\1/' 
fi

