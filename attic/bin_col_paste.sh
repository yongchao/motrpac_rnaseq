#!/bin/bash

#future plan is replace row_paste.awk with the script
set -euo pipefail
optskip=0
optinfo=1
while getopts s:i:c:h o 
do      
    case "$o" in
	s) optskip="$OPTARG";;
	i) optinfo="$OPTARG";;
	c) optcol="$OPTARG";;
	h) echo "Usage: $0 [-h] [-s skip] [-c col] [-i info]"
	   echo '-h: print help'
	   echo '-s: The number of lines to skip'
	   echo '-c: which column to collect, starts with one'
	   echo '-i: which gave the id, or row names'
	   exit 0;;
    esac
done

