#This script needs to be run at the terminal interactively
#It will install python2 and python3 respectively at the folder $conda/python2 and $conda/python3 and other bioninfomratics softwares/packages
#required by motrpac
#The actual python2 and python3 exceutable files will be $conda/python2/bin/python2 $conda/python3/bin/python3

#$conda will be eventually pointed by MOTRPAC_conda
#  The variable conda needs to be defined here, for example
#  export conda=/sc/orga/projects/sealfs01a/conda #this is an example for sinai
cd $conda

#During both installations, always answer no for the question "Do you wish the installer to prepend the Miniconda2 install location ..." (or the
#the question "Do you wish the installer to initialize Miniconda2 in your ..."

#We want to load the conda system on the fly rather than in all bash

#Installing at the folder $conda/python2
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -p python2

#Installing at $conda/python3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p python3

#We will mostly use python3, but sometimes it requires python2, so appending the python3 to the path
export PATH="$conda/python3/bin:$PATH"
#add channels
conda config --prepend channels conda-forge
conda config --prepend channels bioconda


#Uncomment the line below to find conflict pacakges
#conda install python=3.7.1 snakemake=5.3.0 star=2.6.1b cutadapt=1.18 picard=2.18.16 samtools=1.9 multiqc=1.6 rsem=1.3.1 bowtie2=2.3.4.3 fastqc=0.11.8 r-base=3.5.1

#The conflicts tested on Nov 20, 2018
#1. snakemake requires python=3.6.6 so python=3.7.1 can not be used.
#2. rsem relies on samtools 1.3.1 and r-base=3.4.1 so samtools=1.9 and r-base=3.5.1 can not be used
#3. The final version as of Dec 7, 2018

conda install \
      python=3.6.6 \
      snakemake=5.3.0\
      star=2.7.0d\
      cutadapt=1.18 \
      picard=2.18.16 \
      samtools=1.3.1 \
      r-base=3.4.1 \
      rsem=1.3.1 \
      multiqc=1.6 \
      bowtie2=2.3.4.3\
      fastqc=0.11.8 \
      bismark=0.20.0\
      subread=1.6.3\
      ucsc-gtftogenepred=366\
      gawk=4.2.1
