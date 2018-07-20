# RNA-seq pipeline usage
The usage for the RNA-seq pipeline
- running the load_motrapac.sh as specified in the main [README.md](README.md) file
- Each fastq files in the `fastq` subfolder should have `_R1.fastq.gz` or `_R2.fastq.gz` suffix.
- go to the `fastq` folder, run the command prep_config.sh to prepare the file `sample_info.txt` and a preliminary file `config.yaml`
- change the genome folder as neccessary in the file `config.yaml`
- go to top folder fo the project, and the run the command `snakemake -f $MOTRPAC_ROOT/rna-seq.snakemake` for debugging
- If things are working OK, kill the above command and the run the command `Snakemake_lsf.sh -- -f $MOTRPAC_ROOT/rna-seq.snakemake` for production run.

## Python software versions and location
- python/2.7.15 at conda/python2
- python/3.6.3 at conda/python3

## Softwares in conda/python2
- Most packages should be in python3 unless we don't have a choice

## Softwares in conda/python3
- It seems not supporting multi-threads
- rsem/1.2.21 (this version is too low) <1.3.0
- R/3.4.1 (this was installed by rsem, it may be removed later) <3.5.1
- bowtie2/2.2.4 <2.2.8
- multiqc/1.6a0
- cutadapt/1.16

## Software provied by loading minerva modules
- fastqc/0.11.7
- star/2.5.4b
- samtools/1.4.1
- java/1.8.0_111
- picard/2.7.1
