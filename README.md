# snakemake implementation of MoTrPAC RNA-seq pipeline
**Contact:** Yongchao Ge (yongchao.ge@mssm.edu)

* [MoTrPAC RNA-seq MOP (web view version 2.0)](https://docs.google.com/document/d/e/2PACX-1vRFurZraZfxfMd5BWfIQEnETlalDNjQPyMjS7TCTgc3MMlMtB_-tmJfEK7lmRV7GD30I7R9-ISX3kuM/pub)
* This README file and the MOP\_details files were also consulted with Nicole Renee Gay's implementation https://bitbucket.org/nicolerg/motrpac_rna_mop 

# A. External softwares installation and bash environmental setup

## A.1 Conda installation 
We are heavily relying on conda to install/update many bioinformatics software with the same fixed versions. Most updated softwares are available at conda https://conda.io/miniconda.html . 
* Following the instructions at [conda\_install.sh](bin/conda_install.sh) to install the python2 and python3 under the conda root folder `$conda` (a user defined path to install conda). The actual python2 and python3 executable files will be `$conda/python2/bin/python2` and  `$conda/python3/bin/python3`.
* The last command at the file `conda_install.sh` installs the specified versions of software packages and their dependency packages.
  ```bash
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
      fastqc=0.11.8\
      bismark=0.20.0\
      subread=1.6.3\
      ucsc-gtftogenepred=366\
      gawk=4.2.1
  ```

## A.2 Bash environment setup
We rely on the same set-up of conda installation so that the dependency softwares/packages (with the same versions) are portable. We need to export the environmental variables `PATH`, `MOTRPAC_root`,`MOTRPAC_conda`,`MOTRPAC_refdata`
* `MOTRPAC_root` is the root folder of the github code 
* `MOTRPAC_conda` is conda installation folder `$conda` mentioned above 
* `MOTRPAC_refdata` is the genome reference folder
* `export $(bin/load_motrpac.sh)`will export all of the above environment variables for Sinai people
* `export $(bin/load_motrpac.sh -c $conda -r $refdata)` for non-sinai people, where `$conda` is the conda installation folder as mentioned in section A.1 and `$refdata` is the folder that will contain the genome data as mentioned in section A.3 below.
* If the above `export` may not work, try to run the command directly `bin/load_motrpac.sh` to identify errors, and then run the corresponding export commands. 

## A.3 Download the genome source and build the refdata
* Follow the commands in [source\_data.sh](bin/source_data.sh) to download the genome source data (fa and gtf from gencode and ensembl) and also build the bowtie2\_index for the miscellaneous small data (globin and rRNA)
* Running the snakefile [genome\_index.snakefile](genome\_index.snakefile) to build the genome index for each genome folder that was downloaded by [source\_data.sh](bin/source_data.sh). The following is an example hg38\_gencode\_v30
  ```bash
  cd  $MOTRPAC_refdata/hg38_gencode_v30
  snakemake -s $MOTRPAC_root/genome_index.snakefile
  #The above snakemake command can use more than one CPU core to speed it up as below
  #snakemake -j $NUMBEER_OF_CPUS -s $MOTRPAC_root/genome_index.snakefile
  #The computation be further sped up by submitting the jobs to clusters, see section B.2 on the details.
  ```
* The following gives the implementation technical details of the above commands and *these commands shouldn't run separately*.  
  * Human data is on `hg38_gencode_v30` from gencode
    ```bash
	ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.primary_assembly.annotation.gtf.gz
	ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.primary_assembly.genome.fa.gz
    ```
  * Rat data is on `rn6_ensembl_r96` from ensembl
    ```bash
	ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
	ftp://ftp.ensembl.org/pub/release-96/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.96.gtf.gz
    ```
  * The gtf and fa files from ensembl have been modified to have "chr" as part of the chromosome name as in gencode data, see [fixchr4ensembl.sh](bin/fixchr4ensembl.sh) for details
  * The gtf file is sorted, details can be seen in [genome.sh](bin/genome.sh)
  * Different types genome reference index was built with the following command and the details on these commands will be described in [MOP_details.md](MOP_details.md)
	* [star\_index](bin/star_index.sh)
	* [bowtie2\_index](bin/bowtie2_index.sh)
	* [rsem\_index](bin/rsem_index.sh)
	* [qc53\_ref](bin/qc53_ref.sh)
	* [bismark\_index](bin/bismark_index.sh), only used for RRBS data, not used for RNA-seq data.


# B. Running the pipeline

## B.1 bcl2fastq
Details on setting-up the sequencing parameters for NuGEN are described in Section 1 of the MOP.
* [bcl2fastq.sh](bin/bcl2fastq.sh) generates three fastq files, R1, I1 and R2. 
  All of the fastq files can be saved with `${SID}_R1.fastq.gz`, `${SID}_I1.fastq.gz` and `${SID}_R2.fastq.gz` under folder `fastq_raw`, where `SID` is the sample id. 
* Using the mask in bcl2fastq command `--use-bases-mask Y*,I8Y*,I*,Y*` and options `--mask-short-adapter-reads 0 --minimum-trimmed-read-length 0`
* Joining fastq files from all lanes and renaming the files, note that the old `R2` becomes new `I1` and old `R3` becomes new `R2`. 
  ```bash
    cat  $fastq_folder/${SID}_S*_L00?_R1_001.fastq.gz  >${SID}_R1.fastq.gz
    cat  $fastq_folder/${SID}_S*_L00?_R2_001.fastq.gz  >${SID}_I1.fastq.gz
    cat  $fastq_folder/${SID}_S*_L00?_R3_001.fastq.gz  >${SID}_R2.fastq.gz
  ```

## B.2 Run the snakemake program
* In a work folder, a subfolder `fastq_raw` contains the fastq files of all samples `${SID}_R1.fastq.gz`, `${SID}_I1.fastq.gz` and `${SID}_R2.fastq.gz`.
* Make sure the `MOTRPAC_root`, `PATH` and other environmental variables have been setup correctly according to section A.2

### B.2.1 Run the snakemake locally  
* Run the command locally to debug possible problems below for the human genome  
  `snakemake -s $MOTRPAC_root/rna-seq.snakefile`  
* If the data is for rat samples, run the command below for the rat genome `rn6_ensembl_r96`  
  `snakemake -s $MOTRPAC_root/rna-seq.snakefile --config genome=rn6_ensembl_r96`  
  
### B.2.2 Run the snakemake on a cluster
If the snakemake is running OK locally, then submit the snakemake jobs to the cluster. This is only necessary for large jobs. These scripts were written for Sinai LSF and Stanford SLURM job submission systems. Other cluster job submission systems may need to write their own scripts and configuration files (`$MOTRPAC_root/config/`). 

**For the Sinai LSF job submission system:**  
```
Snakemake_lsf.sh -- -s $MOTRPAC_root/rna-seq.snakefile --config genome=rn6_ensembl_r96
```

**For a SLURM job submission system (e.g. SCG Informatics Cluster):**  
* Where `${genome}` is defined as `hg38_gencode_v30` for human RNA-seq data or `rn6_ensembl_r96` for rat RNA-seq data, run the following command to run the pipeline with `sbatch`:  
```
$MOTRPAC_root/bin/Snakemake_slurm.sh ${genome} ${outdir} 
```
* Change `SBATCH` options as needed in `$MOTRPAC_root/config/slurm.json`.   

# *Code implementation philosophy 
* The default uses `python3`, while a few tools relies on `python2`. Both python3 and python2 co-exists peacefully by calling `python2` for python2 specific scripts.
* Each individual component was implemented as an independent bash script file
* The bash script follows the [strict bash mode](http://redsymbol.net/articles/unofficial-bash-strict-mode/) with the setting of `set -euo pipefail`
* The pipeline works for single ends file (with only `_R1` file) or paired-ends file (with `_R1` and `_R2` files) or with UMI (with `_I1` file)
* The default option is for the paired-end files with UMI of MOTRPAC RNA-seq with NuGEN UDI
* The file [MOP\_details.md](MOP_details.md) describes UMI paired ends setting (the same as what the MOP described), but the actual bash script handles more complicated situations for different species and different types of data as mentioned right above.
* With the appropriate setting of the parameters, this pipeline can also work with other RNA-seq data
* Some codes are also part of other MOTRPAC projects (for example, RRBS-seq data) for building the genome references.

# *The implementation details
* The implementation details on the RNA-seq MOP can be seen in [MOP\_details.md](MOP_details.md). 
* The commands in [MOP\_details.md](MOP_details.md) don't need to be run separately as everything has already been taken care of by the snakemake commands in the above Section B.2 

