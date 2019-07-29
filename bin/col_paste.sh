#!/bin/bash 
set -euo pipefail
#col_paste.sh -s 1 -c 5 -i 1 sample_info.txt >out.txt
#to replace row_paste.awk with the script
#opefully this is going to be faster and more streamlineds

##samples.txt is tab delimited file with at least two columns, column 1 is the the file name and column 2 is the sample.id
#if it has only one column, then the file name is also the sample.id
##collect all of the data from many samples
##Assuming one line for the file header
##No space is allowed in file sample_info.txt, but it is OK for otehr entries
##And Tab delimited file

mkdir -p tmpdir
tmpdir=$(mktemp -d -p $MOTRPAC_tmp col_paste.XXXX) #create a temp folder to be easier to be removed
optskip=0
optinfo=1
optcol=0
while getopts s:i:c:h o 
do      
    case "$o" in
        s) optskip="$OPTARG";;
        i) optinfo="$OPTARG";;
        c) optcol="$OPTARG";;
        h) echo "Usage: $0 [-h] [-s skip] [-c col] [-i info]"
           echo '-h: print help'
           echo '-s: The number of lines to skip (s==-1 for no file header), default 0'
           echo '-c: which column to collect, starts with one, default the last column'
           echo '-i: which gave the id, or row names, default 1'
           exit 0;;
    esac
done
shift $((OPTIND-1))


infile=$1
files=(`cut -f 1 $infile`)
NFinfile=$(awk 'NR==1{print NF}' $infile)
if ((NFinfile >= 2)); then
    samples=(`cut -f 2 $infile`)
else
    samples=("${files[@]}") 
fi

tmpfiles=""
rowid0=$tmpdir/rowid0 #the row id

for i in "${!samples[@]}" #zero based index for bash array
do
    
    sid=${samples[$i]}
    file=${files[$i]}
    rowid=$tmpdir/rowid
    awk -F $'\t' 'NR>'$((optskip+1))'{print $'$optinfo'}' $file  >$rowid 
    if ((i == 0)); then
        NF=$(awk -F $'\t' 'NR=='$((optskip+2))'{print NF}' $file)
	if ((optcol==0)); then ((optcol=NF));fi
	if ((optskip>=0)); then
	    printf $(awk -F $'\t' 'NR=='$((optskip+1))'{print $'$optinfo'}' $file)
	else
	    printf "field"
	fi
	tmpfiles=$rowid
	cp $rowid $rowid0
    else
	cmp -s $rowid $rowid0 || {
	    echo "The rowids do not match for sample $sid\n" >& 2
	    exit 1
	}
    fi
    #for all samples
    printf "\t"$sid
    tmpfile=$tmpdir/$sid
    tmpfiles="$tmpfiles $tmpfile"
    awk -F $'\t' 'NR>'$((optskip+1))'{print $'$optcol'}' $file  >$tmpfile
done
printf "\n"
paste $tmpfiles
rm -rf $tmpdir
