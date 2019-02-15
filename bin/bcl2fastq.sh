#This script is best to run interactively to find where the problem is.

#Please note that SampleSheet.csv needs to have the adpater setting removed
#as we are going to use cutadapt to remove the adapter anyway.

mkdir bcl2fastq
cd bcl2fastq

# $bcl_folder is defined to be the folder that contains the raw data of the bcl run folder
# It can be defined here by
# export bcl_folder=/path/to/bcl_raw_data

bcl2fastq --sample-sheet SampleSheet.csv \
	  --use-bases-mask Y*,I8Y*,I*,Y* \
	  -R  $bcl_folder \
	  --mask-short-adapter-reads 0 \
	  --minimum-trimmed-read-length 0 \
	  -o . 
#Please note that the options are necessary --mask-short-adapter-reads 0, --minimum-trimmed-read-length 0
#for the index reads, otherwise all index reads are masked	  

#The following command joins the fastq files from bcl2sfastq to fastq_raw folders
#rename R2 to I1 and R3 to R2.

#getting the samples from sample sheet

samples=$(sed -n '/^\[Data\]/,$p' SampleSheet.csv |tail -n +3|cut -f 2 -d ,)
cd ..
mkdir -p fastq_raw
cd fastq_raw
for SID in $samples; do
    fastq_folder=$(dirname $(find ../bcl2fastq -name "${SID}_S*_L00?_R1_001.fastq.gz" |head -1))
    cat  $fastq_folder/${SID}_S*_L00?_R1_001.fastq.gz  >${SID}_R1.fastq.gz
    cat  $fastq_folder/${SID}_S*_L00?_R2_001.fastq.gz  >${SID}_I1.fastq.gz
    cat  $fastq_folder/${SID}_S*_L00?_R3_001.fastq.gz  >${SID}_R2.fastq.gz
done

