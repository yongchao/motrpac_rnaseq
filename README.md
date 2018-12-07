# snakemake implementation of MoTrPAC RNA-seq pipeline
**Contact:** Yongchao Ge (yongchao.ge@mssm.edu)

* The initial MoTrPAC RNA-seq MOP: https://docs.google.com/document/d/1oz8jZAY9Rq4uqenp-0RMkQBhMtjsMLKlxfjlmRExSQ0/edit?ts=5b04a52e#
* This README file was also consulted with Nicole Renee Gay's implementation https://bitbucket.org/nicolerg/motrpac_rna_mop 

# A. External softwares installation and bash environmental setup

## A.1 Conda installation 
We are heavily relying on conda to install/update many bioinformatics software with the same fixed versions. Most updated softwares are available at conda https://conda.io/miniconda.html . 
* Install the python2 and python3 under the conda root folder `$conda` following the instructions at [conda\_install.sh](bin/conda_install.sh)
* The last command at the file `conda_install.sh` installs the specified versions of software packages and their dependency packages.
  ```bash
  conda install \
      python=3.6.6 \
      snakemake=5.3.0\
      star=2.6.1b\
      cutadapt=1.18 \
      picard=2.18.16 \
      samtools=1.3.1 \
      r-base=3.4.1 \
      rsem=1.3.1 \
      multiqc=1.6 \
      bowtie2=2.3.4.3\
      fastqc=0.11.8\
      bismark=0.20.0\
      subread=1.6.2
  ```

## A.2 Bash environments setup
We rely on the same set-up of conda installation so that the dependency softwares/packages (with the same versions) are portable. We need to export the environmental variables `PATH`, `MOTRPAC_root`,`MOTRPAC_conda`,`MOTRPAC_refdata`
* `MOTRPAC_root` is the root folder of the github code 
* `MOTRPAC_conda` is conda installation folder `$conda` mentioned above
* `MOTRPAC_refdata` is the genome data index folder
* `export $(bin/load_motrpac.sh)`will export the all of the above environment variables for Sinai people
* `export $(bin/load_motrpac.sh -c $conda -r $refdata)` for non-sinai people. 

## A.3 Download the genome source and build the refdata
* Download the genome source data (fa and gtf from gencode and ensembl) and also build the bowtie2\_index for the miscellaneous small data (globin and rRNA)
  [source\_data.sh](bin/source_data.sh)
  * Human data is on `hg38_gencode_v29` from gencode
    ```bash
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
    ```
  * Rat data is on `rn6_ensembl_r94` from ensembl
    ```bash
    ftp://ftp.ensembl.org/pub/release-94/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gtf.gz
    ftp://ftp.ensembl.org/pub/release-94/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
    ```
* Make sure the gtf file and fa file from ensembl have "chr" as part of the chromosome name as in gencode data (see [fixchr4ensembl.sh](bin/fixchr4ensembl.sh))
* Sort the gtf file accordingly, see [genome.sh](bin/genome.sh)
* Build the genome reference data, the details on these commands will be described latter.
  * [star\_index](bin/star\_index.sh)
  * [bowtie2\_index](bin/bowite2\_index.sh)
  * [rsem\_index](bin/rsem\_index.sh)
  * [qc53\_ref](bin/qc53\_ref.sh)
* For each genome folder that was downloaded by [source\_data.sh](bin/source\_data.sh). The genome data and index for the refdata can be built at the refdata with the command, for example hg38\_gencode\_v29
  ```bash
  cd  $MOTRPAC_refdata/hg38_gencode_v29
  snakemake -s $MOTRPAC_root/genome_index.snakefile
  ```

# B. Running the pipeline

## B.1 bcl2fastq
* [bcl2fastq.sh](bin/bcl2fastq.sh) generates three fastq files, R1, I1 and R2.
  All of the fastq files can be softlinked with `${SID}_R1.fastq.gz`, `${SID}_I1.fastq.gz` and `${SID}_R2.fastq.gz`, where `SID` is the sample id. And all of these files are saved in the sub folder `fastq_raw`
* Using the mask in bcl2fastq command `--use-bases-mask Y*,I8Y*,I*,Y*`
* Renaming the three files `_R1`, `_R2`, `_R3` to `_R1`, `_I1` and `_R2` and getting rid off the non-meaningful file names. 
  ```bash
  find . -name "*_R2_001.fastq.gz" |
  awk '{print "mv "$0" "gensub("_R2_001.fastq.gz$","_I1_001.fastq.gz",1)}' |bash -x
  find . -name "*_R3_001.fastq.gz" |
  awk '{print "mv "$0" "gensub("_R3_001.fastq.gz$","_R2_001.fastq.gz",1)}' |bash -x
  ```
* Getting rid off the non-meaning names.
  ```bash
  find ../bcl2fastq -name "*_001.fastq.gz" |
  awk 'BEGIN{FS="/"};{print "ln -s "$0" "gensub(/_S[0-9]+_L00[1-8](_[IR][12])_001/,"\\1",1,$(NF))}' |bash -x
  ```

## B.2 Run the snakemake program
* In a work folder, a subfolder `fastq_raw` contains the fastq files of all samples `${SID}_R1.fastq.gz`, `${SID}_I1.fastq.gz` and `${SID}_R2.fastq.gz`.
* Make sure the `MOTRPAC_ROOT` and `PATH` have been setup correctly according to section A.2
* Run the command locally to debug possible problems below for the human genome  
  `snakemake -s $MOTRPAC_root/rna-seq.snakefile`
* If the data is for rat samples, run the command below for the rat genome `rn6_ensembl_r94`  
  `snakemake -s $MOTRPAC_root/rna-seq.snakefile --config genome=rn6_ensembl_r94`
* If the snakemake is running OK locally, then submit the snakemake jobs to the cluster. This is only necessary for large jobs. This script was writte for Sinai LSF jobs submission system. Other cluster job submission system may need to write their own script.  
  `Snakemake.lsf -- -s $MOTRPAC/rna-seq.snakefile --config genome=$genome`

## B.3 Code implementation philosophy 
* The default is using `python3`, while a few tools relies `python2`. Both python3 and python2 co-exists peacefully by calling `python2` for python2 specific tools.
* Each individual component was implemented as an independent bash script file
* The bash script follows the [strict bash mode](http://redsymbol.net/articles/unofficial-bash-strict-mode/) with the setting of `set -euo pipefail`
* The pipeline works for single ends file (with only `_R1` file) or paired-ends file (with `_R1` and `_R2` files) or with UMI (with `_I1` file)
* The default option is for the paired-end files with UMI of MOTRPAC RNA-seq with NuGEN UDI
* The implementation details in this document only descrirbes the UMI paired ends setting, but the actual bash script handles more complicates siutations for different species and different types of data as mentioned right above.
* With the appropriate setting of the parameters, this pipeline can also work with other RNA-seq data
* Some codes are also part of other MOTRPAC projects (for example, RRBS-seq data)

# C. Pre-alignment sample processing and QC

## C.1 Run fastqc on the fastq file in fastq_raw folder
* The fastq command details in [fastqc.sh](bin/fastqc.sh)
  ```bash
  fastqc -o $odir $fqfile 
  ```
* `$odir`:  `fastqc_raw`
* `$fqfile`: one of the fastq files in `${SID}_R1.fastq.gz` and `${SID}_R2.fastq.gz` under folder `fastq_raw`

## C.3 Attach UMI from index file to read files
In order to allow the the later data to track the UMI for each read, one stragety is to attach the UMI form the index fastq file `_I1` to `_R1` and `_R2`. The implementation is in [UMI_attach.awk](bin/UMI_attach.awk). We will apply this awk command for all of `_R1` and `_R2` files in the folder `fastq_raw` and save them into folder `fastq_attach`

## C.3 Adapter trimming
For each paired-end fastq files (`${SID}_R1.fastq.gz` and `${SID}_R2.fastq.gz`) under folder `fastq_attach`, remove adapters using cutadapt. Eliminate reads that are too short after removing adapters and save trimmed fastq files into folder `fastq_trim`. Save reads that were trimmed because they were too short into the folder `fastq_trim/too_short`. The details are in [trim.sh](bin/trim.sh). The following will be for the paired input fastq files. The implementation also considers single ended fastq file.
```bash
#Following trim_galore approach, taking the first 13 bases
#This is due to low quality at 3 prime end and also NuGEN may not take the complete universtal adapter
index_adapter=AGATCGGAAGAGC
univ_adapter=AGATCGGAAGAGC
R1=${SID}_R1.fastq.gz
R2=${SID}_R2.fastq.gz
cutadapt \
    -a $index_adapter \
    -A $univ_adapter \
    -o fastq_trim/$R1 \
	-p fastq_trim/$R2\
    -m 20 \
    --too-short-output fastq_trim/tooshort/$R1\
    --too-short-paired-output fastq_trim/tooshort/$R2\
	fastq_raw/$R1 fastq_raw/$R2
```
## C.4  Run fastqc on all of the fastq files in fastq\_trim folder
* In the command [fastqc.sh](bin/fastqc.sh), specifying `$odir` to be `fastqc` and `$fqfile` to be one of the fastq files in `${SID}_R1.fastq.gz` and `${SID}_R2.fastq.gz` under the folder `fastq_trim`

## C.5 MultiQC for pre-aligned data, collecting the QC metrics on the fastqc, fastqc\_raw and fastq\_trim folders
The implementation is in the rule `pre_align_QC` of the `rna-seq.snakefile` and [multiqc.sh](bin/multiqc.sh)
```bash
multiqc\
	-d
	-f \
	-n pre_align \
	-o multiqc \
	fastqc_raw fastq_trim fastqc
```

# D. Alignments 
Please note all D.[0-9]A parts are to describe the implmentation details for generating the reference. All of these references can be generated in one single command as in the last item of section A.3 and only needs to be done once.
## D.1A Generate the star index
The star index is implemented in [star\_index.sh](bin/star\_index.sh). `genome.fa` and `genome.gtf` are the cleaned-up fastq and gtf files of human or rat genome (the clea-up was decsried in section A.3). Note that we use the default value 100, which works as good as sequence length specific star index. 
```bash
STAR --runThreadN $threads\
     --runMode genomeGenerate\
     --genomeDir star_index\
     --genomeFastaFiles genome.fa\
     --sjdbGTFfile genome.gtf\
     --sjdbOverhang 100\
     --outFileNamePrefix star_index/
```
## D.1B Align paired fastq files in the fastq\_trim folder and save the results in star\_align folder
For each paired fastq files `${SID}_R1.fastq.gz` and `${SID}_R2.fastq.gz` under folder `fastq_trim`, align the reads to the genome ([star\_align.sh](bin/sta\r_align.sh)). 
Parameters:
* `gdir`: the folder that contains the file `genome.fa` and `genome.gtf` and `star_index` and `rsem_index` and `bowtie_index`. For human, it will be `$MOTRPAC_refdata/hg38_gencode_v29`
```bash
mkdir -p star_align/$SID
STAR  --genomeDir $gdir/star_index\
      --sjdbOverhang  99\
      --readFilesIn fastq_trim/${SID}_R1.fastq.gz fastq_trim/${SID}_R2.fastq.gz\
      --outFileNamePrefix star_align/${SID}/\
      --readFilesCommand zcat \
      --outSAMattributes NH HI AS NM MD nM\
      --outFilterType BySJout\
      --runThreadN $threads\
      --outSAMtype BAM SortedByCoordinate\
      --quantMode TranscriptomeSAM
```
## D.1C Quantify by featureCounts
The implementation details are in the file [featureCounts.sh](bin/featureCounts.sh). The essential part is
```bash
featureCounts -a $gdir/genome.gtf -o featureCounts/$SID -p -M --fraction star_align/${SID}/Aligned.sortedByCoord.out.bam
```

## D.2A Prepare RSEM reference
It needs to be done once [rsem\_index.sh](bin/rsem\_index.sh). 
```bash
rsem-prepare-reference --gtf genome.gtf genome.fa rsem_index/genome
```
## D.2B RSEM quantification 
Run quantification with the sorted `star_algin/$SID/Aligned.toTranscriptome.sorted.bam` files. For paired fastq files, the implementation from [rsem.sh](bin/rsemh.sh) is like below
```bash
rsem-calculate-expression \
    --paired-end \
    -p $threads\
    --no-bam-output\
    --forward-prob 0.5\
    --seed 12345\
    --bam star_align/$SID/Aligned.toTranscriptome.out.bam\
    $gdir/rsem_index/genome\
    rsem/$SID
```
## D.3A Prepare bowtie2 index for globin and rRNA and ERCC and phix
This is only needed to be done once and is implemented in the file [source\_data.sh](bin/source\_data.sh). The following is the essential part for human rRNA bowtie2 index
```
bowtie2-build hg_rRNA.fa $hg_rRNAref
```
## D.3B Quantify the alignment percentage for globin and rRNA and ERCC and phix
The implemention is at [bowtie2.sh](bin/bowtie2.sh). This is how the code will work for the human rRNA. 
```bash
bowtie2 -p $threads \
	-1 fastq_trim/{SID}_R1.fastq.gz \
	-2 fastq_trim/{SID}_R2.fastq.gz \
	-x $hg_rRNAref 
	--local \
	-S $sam >rRNA/{SID}.txt
```
### D.4A prepare refFlat file for Picard CollectRnaSeqMetrics
* refFlat files are required for  Picard CollectRnaSeqMetrics. The refFlat file is generated from the GTF file of each genome folder under `MOTRPAC\_refdata`. The implementation details can be seen in [qc53\_ref.sh](bin/q53\_ref.sh)
### D.4B Collect RNA-seq metrics with Picard CollectRnaSeqMetrics
Compute post alignment QC metrics, including % mapped to coding, intron, intergenic, UTR, and % correct strand, and 5’ to 3’ bias.  
The command for computing the CollectRnaSeqMetrics is in [qc53\_ref.sh](bin/q53\_ref.sh) and the essential part is
```bash
picard CollectRnaSeqMetrics\
     I=$bam \
     O=qc53/${SID}.RNA_Metrics\
     MINIMUM_LENGTH=50 \
     REF_FLAT=$gdir/qc53_ref/genome_flat.txt\
     STRAND_SPECIFICITY=$strand \
     RRNA_FRAGMENT_PERCENTAGE=0.3 \
     >& qc53/log/${SID}.log
```
## D.5 Using the UMI design to remove PCR duplicates
We will use the UMI and the aligned information to find out the duplicates and remove the duplicates. The details are in [UMI.sh](bin/UMI.sh). It involves two steps
* Use the samtools to find the uniquely mapped reads from the bam file, see file [bam\_uniq.sh](bin/bam_uniq.sh)
* Use the nugen python2 script file `nudup.py` to compute the duplicattion rate and the duplicates removed bam file
* This process is applied for the all bam files in the folder `star_align/$SID` and save the results to `star_align/UMI_$SID`
* For the bam files in `star_align/UMI_$SID`, we will apply the computation D.2B and D.4B

## D.6 MultiQC for post-aligned data, collecting the QC metrics on the star_align rsem and featureCounts folders
The implementaiton is in the rule `pre_align_QC` of the `rna-seq.snakefile` and [multiqc.sh](bin/multiqc.sh)
```bash
multiqc \
	-d \
	-f \
	-n pre_align \
	-o multiqc \
	fastqc fastqc fastqc_raw
```
# E Compile important metrics from the MultiQC output and other log files
The R script [qc.R](qc.R) collects all of the important the QC metrics from multiQC output and other log files. All of the metrics have been saved
into the file `qc_info.txt` after the pipeline finishes.

## E.1 fastq metrics (raw and trimmed), collected from pre-alignment QC (see section C.5).
For paired ends fastq files, the metrics has been averaged between the two files to get a single metric.
* reads_raw: number of reads in raw fastq files (in folder `fastq_raw`)
* %trimmed: % trimmed reads 
* %trimmed_bases: % trimmed bases
* reads: number of reads in the main fastq files (in folder `fastq_trim` when trimmed or in folder `fastq` when no trimming)
* %GC: percentage of the GC content
* %dup_sequence: percentage of duplicated sequences

## E.2 Contamination and presence of other common spike-ins in the trimmed fastq files (see Section D.3B)
* %rRNA
* %globin
* %ERCC
* %phix
## E.3 The PCR duplication rate marked by PICARD MarkDuplicates
The implementation is in the file [mark\_dup.sh](bin/mark_dup.sh)
* %dup_picard
## E.4 The PCR duplication rate assessed with UMI technique (see section D.5)
* %umi_dup

## E.5 Post\-alignment RNA-seq metrics from Picard's CollectRnaSeqMetrics output files, also for UMI duplicated removed files
* The important metrics obtained from the Picard's CollectRnaSeqMetrics output file. 
  * %coding: % of bases mapped to coding
  * %utr: % of bases mapped to UTR
  * %intronic: % of bases mapped to intronic
  * %intergenic: % of bases mapped to intergenic
  * %mrna: % of bases mapped to mRNA
  * median\_5'\_3'\_bias:	median 5prime to 3 prime bias
* Similar metrics obtained for the duplicated removed reads based in UMI  
  * %coding
  * %umi_utr
  * %umi_intronic
  * %umi_intergenic
  * %umi_mrna
  * umi_median\_5'\_3'\_bias

## E.6 Alignment metrics from STAR's ${SAMPLE}\_Log.final.out files
* avg\_input\_read\_length
* uniquely\_mapped
* %uniquely\_mapped	
* avg\_mapped\_read\_length	
* num\_splices	
* num\_annotated_splices	
* num\_GTAG\_splices	
* num\_GCAG\_splices	
* num\_ATAC\_splices	
* num\_noncanonical\_splices	
* %multimapped	
* %multimapped\_toomany	
* %unmapped\_mismatches	
* %unmapped\_tooshort	
* %unmapped\_other	
* %chimeric
## E.7 Alignment metrics from RSEM ${SAMPLE}.stat/${SAMPLE}.cnt files, also for UMI duplicated removed files
*Is this section necessary?*
* N\_ALIGNABLE: 
* N\_UNIQUE: 
* N\_MULTI: 
* N\_UNCERTAIN:
* N\_TOTAL\_ALIGNMENTS: 

## E.7 Alignment metrics,also for UMI duplicated removed files
See the file [bam_chrinfo.sh](bin/chrinfo.sh) for details
* for original bam file, `star_align/${SID}`. In the calculation, we only considered the primary alignment.
  * %chrX:  % of reads mapped to ChrX
  * %chrY:  % of reads mapped to ChrY
  * %chrM:  % of reads mapped to ChrM
  * %chrAuto:  % of reads mapped to autosomes
  * %contig:  % of reads mapped to contigs
* similiar metrics for UMI duplicated removed files. The bam files from `star_align/UMI_${SID}`
  * %umi\_chrX:
  * %umi\_chrY:  
  * %umi\_chrM: 
  * %umi\_chrAuto: 
  * %umi\_contig:

# F. Flag problematic samples

*These metrics are very liberal to remove samples of very low quality. Additional samples might be removed downstream using more stringent criteria.*

Flag samples with 

1. Low RIN scores, e.g. RIN < 6
2. Low number of sequenced reads after trimming adaptors, e.g. samples with < 20M reads
3. Abnormal GC content, e.g. GC% > 80% or < 20%
4. High percentage of rRNA, e.g. rRNA > 20%
5. Low number of mapped reads, e.g. < 60%
6. Low number of % exonic read, % exonic read < 50%
7. Mapped read count < 50% of average mapped read count per sample  

# G. Post-Quantification QC 

*These steps should be updated when new samples become available because they are depend on all samples.*

1. Confirm sample concordance with other omics genotypes from same individual using VerifyBAMID.
2. Obtain filtered and normalized expression data and filter outliers:  

  * Remove genes with 0 counts in every sample. 
  * Filter lowly expressed genes (for QC purposes only, definition of expressed gene might change for each analyses): genes with mean counts < 5 and zero counts in more than 20% of samples.
  * Perform library size correction and variance stabilization using DESeq2 or other methods.
  * Compute D-statistic for each sample based on filtered and normalized expression data, i.e. average Spearman’s correlation with other samples from the same tissue and time point. Flag samples with D<.80.
  * Perform PCA based on filtered and normalized expression data across and within tissues. Mark outlier samples (outside +/- 3 sd).  
