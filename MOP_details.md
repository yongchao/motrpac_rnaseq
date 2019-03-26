# snakemake implementation of MoTrPAC RNA-seq pipeline
**Contact:** Yongchao Ge (yongchao.ge@mssm.edu)

* The MoTrPAC RNA-seq MOP: https://docs.google.com/document/d/1oz8jZAY9Rq4uqenp-0RMkQBhMtjsMLKlxfjlmRExSQ0/edit?ts=5b04a52e#
* The other file [README.md](README.md) describes on how to run the pipeline
* This document describes the implementation details for the MOP. No commands in this document should run separately.
* The first section in this file starts with section C while sections A and B are in the other file  [README.md](README.md)

# C. Pre-alignment sample processing and QC

## C.1 Run fastqc on the fastq file in fastq_raw folder
* The fastq command details in [fastqc.sh](bin/fastqc.sh)
  ```bash
  fastqc -o $odir $fqfile 
  ```
* `$odir`:  `fastqc_raw`
* `$fqfile`: one of the fastq files in `${SID}_R1.fastq.gz` and `${SID}_R2.fastq.gz` under folder `fastq_raw`

## C.2 Attach UMI from index file to read files
In order to allow the the later data to track the UMI for each read, one strategy is to attach the UMI from the index fastq file `_I1` to `_R1` and `_R2`. The implementation is in [UMI_attach.awk](bin/UMI_attach.awk). We will apply this awk command for all of `_R1` and `_R2` files in the folder `fastq_raw` and save them into folder `fastq_attach`

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
    -p fastq_trim/$R2 \
    -m 20 \
    --too-short-output fastq_trim/tooshort/$R1 \
    --too-short-paired-output fastq_trim/tooshort/$R2 \
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
Please note all D.1A, D.2A, D.3A and D.4A parts are to describe the implementation details for generating the reference. All of these references can be generated in one single command as in the second item of section A.3 of [README.md](README.md) and only needs to be done once.
## D.1A Generate the star index
The star index is implemented in [star\_index.sh](bin/star\_index.sh). `genome.fa` and `genome.gtf` are the cleaned-up fastq and gtf files of human or rat genome (the clea-up was decsried in [README.md](README.md)). Note that we use the default value 100, which works as good as sequence length specific star index. 
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
      --sjdbOverhang  100\
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
It needs to be done once [rsem\_index.sh](bin/rsem_index.sh). 
```bash
rsem-prepare-reference --gtf genome.gtf genome.fa rsem_index/genome
```

## D.2B RSEM quantification 
Run quantification with the sorted `star_algin/$SID/Aligned.toTranscriptome.sorted.bam` files. For paired fastq files, the implementation from [rsem.sh](bin/rsem.sh) is like below
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

## D.3A Prepare bowtie2 index for globin and rRNA and phix
This is only needed to be done once and is implemented in the file [source\_data.sh](bin/source_data.sh). The following is the essential part for human rRNA bowtie2 index
```
bowtie2-build hg_rRNA.fa $hg_rRNAref
```

## D.3B Quantify the alignment percentage for globin and rRNA and phix
The implementation is at [bowtie2.sh](bin/bowtie2.sh). This is how the code will work for the human rRNA. 
```bash
bowtie2 -p $threads \
	-1 fastq_trim/{SID}_R1.fastq.gz \
	-2 fastq_trim/{SID}_R2.fastq.gz \
	-x $hg_rRNAref 
	--local \
	-S $sam >rRNA/{SID}.txt
```

## D.4A prepare refFlat file for Picard CollectRnaSeqMetrics
* refFlat files are required for  Picard CollectRnaSeqMetrics. The refFlat file is generated from the GTF file of each genome folder under `MOTRPAC_refdata`. The implementation details can be seen in [qc53\_ref.sh](bin/qc53_ref.sh)

## D.4B Collect RNA-seq metrics with Picard CollectRnaSeqMetrics
Compute post alignment QC metrics, including % mapped to coding, intron, intergenic, UTR, and % correct strand, and 5’ to 3’ bias.  
The command for computing the CollectRnaSeqMetrics is in [qc53\_ref.sh](bin/qc53_ref.sh) and the essential part is
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

## D.5 Using the UMI design to compute the PCR duplicates
* Use the nugen python2 script file `nudup.py` to compute the duplication rate and the duplicates removed bam file, details are seen in [UMI\_dup.sh](bin/UMI_dup.sh)
```bash
	python2  nudup.py \
        -2  \
       -s 8 \
	   -l 8 \
	   --rmdup-only \
	   -o <output>  \
	   -T $tmpdir \
	   <Aligned.sortedByCoord.out.bam>
```

## D.6 MultiQC for post-aligned data, collecting the QC metrics on the star_align rsem and featureCounts folders
The implementaiton is in the rule `pre_align_QC` of the `rna-seq.snakefile` and [multiqc.sh](bin/multiqc.sh)
```bash
multiqc \
    -d \
    -f \
    -n post_align \
    -o multiqc \
	post_align star_align featureCounts rsem 
```
# E Compile important metrics from the MultiQC output and other log files
The R script [qc.R](bin/qc.R) collects all of the important the QC metrics from multiQC output and other log files. All of the metrics have been saved
into the file `qc_info.csv` after the pipeline finishes.

## E.1 fastq metrics (raw and trimmed), collected from pre-alignment (see section C).
For paired ends fastq files, the metrics has been averaged between the two files to get a single metric. Most qc metrics are extracted directly from Section 3.5 unless indicated otherwise.
* reads_raw: number of reads in raw fastq files (in folder `fastq_raw`)
* %adapter_detected: % of reads with adapters detected (extracted from the cutadapt log file)
* %trimmed: % trimmed reads 
* %trimmed_bases: % trimmed bases
* reads: number of reads in the main fastq files (in folder `fastq_trim` when trimmed or in folder `fastq` when no trimming)
* %GC: percentage of the GC content in trimmed fastq files
* %dup_sequence: percentage of duplicated sequences in trimmed fastq files

## E.2 Contamination and presence of other common spike-ins in the trimmed fastq files (see Section D.3B)
* %rRNA
* %globin
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

## E.6 Alignment metrics from STAR's star\_align/${SID}_Log.final.out files
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

## E.7 Alignment metrics
See the file [bam\_chrinfo.sh](bin/bam_chrinfo.sh) for details
* for original bam file, `star_align/${SID}`. In the calculation, we only considered the primary alignment.
  * %chrX:  % of reads mapped to ChrX
  * %chrY:  % of reads mapped to ChrY
  * %chrM:  % of reads mapped to ChrM
  * %chrAuto:  % of reads mapped to autosomes
  * %contig:  % of reads mapped to contigs

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
