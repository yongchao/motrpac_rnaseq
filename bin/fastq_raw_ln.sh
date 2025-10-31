#!/bin/bash -x

# #1 folder
# #2 lanes 12 for lanes 1 and 2, default for all lanes, if 0, not using lane

##For more lane splitting, see the file ~/bash/fastq_lane.sh

set -eu -o pipefail

fd=../bcl2fastq/
if (( $# >0 )); then
    fd=$1
    if [[ ${fd: -1} != "/" ]];then
	fd=$fd/
    fi
fi

##adding the trail on $fd, future plan

if (( $#==2 )); then
    if (($2==0)); then
	lanes=""
	L=""
    else
	lanes=[$2]
	L=_L00[$lanes]
    fi
else
    lanes=[1-8]
    L=_L00$lanes
fi

#chmod -R a-w  ../bcl2fastq

if [[ -e $fd/SampleSheet.csv ]];then
   sed -n '/^\[Data\]/,$p' $fd/SampleSheet.csv |tail -n +3 >tmpfile
   if (( $#==2 )) && [[ $lanes != "" ]]; then
       samples=$(awk -F , '$2~/^'$lanes'/ {print $1}' tmpfile)
   else
       samples=$(cut -f 1 -d , tmpfile)
   fi
else
    find $fd -name "*${L}_R1_001.fastq.gz">tmpfile
    samples=$(awk 'BEGIN{FS="/"};{sub("_S[0-9]+'$L'_R1_001.fastq.gz$","",$NF); print $NF}' tmpfile)
fi
rm tmpfile

sample1=$(echo $samples| awk 'NR==1{print $1}')

set +e
R3=$(find $fd -name "${sample1}_S*${L}_R3_001.fastq.gz" |grep -P "${sample1}_S[0-9]+${L}_R3_001.fastq.gz"|wc -l)
set -e
#As find can not deal with [0-9]+, has to use grep -P to be more precise

echo me
#R3==0 means no UMI

#currently, it can process the data with or without UMI
for SID in $samples; do
    echo $SID
    fastq_folder=$(dirname $(find $fd -name "${SID}_S*${L}_R1_001.fastq.gz" |head -1))
    ln -s  $fastq_folder/${SID}_S*${L}_R1_001.fastq.gz  ${SID}_R1.fastq.gz
    if [[ $R3 == 0 ]]; then
	ln -s  $fastq_folder/${SID}_S*${L}_R2_001.fastq.gz  ${SID}_R2.fastq.gz
    else
	ln -s  $fastq_folder/${SID}_S*${L}_R2_001.fastq.gz  ${SID}_I1.fastq.gz
	ln -s  $fastq_folder/${SID}_S*${L}_R3_001.fastq.gz  ${SID}_R2.fastq.gz
    fi
done


