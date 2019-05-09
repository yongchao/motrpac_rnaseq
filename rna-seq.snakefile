#Author: Yongchao Ge

#Jan 22, 2018: the version is in sync with the MOP
#Nov 20, 2018: The intial version that is shared with Stanford
#July 20,2018: The initial push to sealfonlab

#The usage are in the file README.md

#Folder structures at root of working folder
#fastq_raw: raw fastq files, with no adpaters removed, probably a softlink to the fastq files in the output folder of bcl2fastq

#For each sample, we have three possible fastq files
#${sid}_R1.fastq.gz, required for all data. 
#${sid}_R2.fastq.gz, required for pairied end reads
#${sid}_I1.fastq.gz, required for NuGEN with UMI for UMI processsing

#configure genome by using --config genome=hg38_gencode_v30 etc
#Updating the link for rn6 and hg38 in $MOTRPAC_refdata so that we don't have to specify the whole name

include: "sample_sub.snakefile"
ruleorder: trim > trim_single

localrules: all,star_align_all,featureCounts_all,rsem_all,samples,fastqc_all

rule all:
    input:
        "log/OK.pre_align_QC",
        "log/OK.post_align_QC",
        "log/OK.qc_all",
        "log/OK.star_align",
        "log/OK.featureCounts",
        "log/OK.rsem"

##As trim's output are dependent on the inout,
rule trim_single:
    input:
        trim_input+"{sample}_R1.fastq.gz"
    output:
        "fastq_trim/{sample}_R1.fastq.gz",
        ok="fastq_trim/log/OK.{sample}"
    log:
        "fastq_trim/log/log.{sample}"
    shell:
        '''
        trim.sh {input} >&{log}
        echo OK>{output.ok}
        '''
rule trim:
    input:
        expand(trim_input+"{{sample}}_{R}.fastq.gz",R=["R1","R2"])
    output:
        expand("fastq_trim/{{sample}}_{R}.fastq.gz",R=["R1","R2"]),
        ok="fastq_trim/log/OK.{sample}"
    priority:
        10 #This is preferred than trim_single, with default priority value of zero
    log:
        "fastq_trim/log/log.{sample}"
    shell:
        '''
        trim.sh {input} >&{log}
        echo OK>{output.ok}
        '''

rule star_align:
    input:
        fastq_info 
    output:
        bam="star_align/{sample}/Aligned.sortedByCoord.out.bam",
        rsem_bam=temp("star_align/{sample}/Aligned.toTranscriptome.out.bam"),
        qc_log="star_align/{sample}/Log.final.out"
    log:
        "star_align/log/{sample}.log"
    threads: 6
    shell:
        '''
        star_align.sh {gdir} {threads} {tmpdir} {input} >&{log}
        '''

rule chr_info:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "star_align/{sample}/chr_info.txt"
    shell:
        '''
        bam_chrinfo.sh {input}

        '''
rule UMI_dup:
    input:"star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:"star_align/{sample}/{sample}_dup_log.txt"
    log:"star_align/{sample}/UMI_dup.log"
    params:
        R2_info
    shell:
        '''
        UMI_dup.sh {input} {params} >&{log}
        '''
rule rsem:
    input:
        "star_align/{sample}/Aligned.toTranscriptome.out.bam"
    output:
        "rsem/log/OK.{sample}"
    log:
        "rsem/log/{sample}.log"
    params:
        R2_info
    threads: 6
    shell:
        '''
        rsem.sh {input} {gdir} {threads} {params} {tmpdir} >&{log}
        echo "Finished rsem" > {output}
        '''
        
rule featureCounts:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "featureCounts/{sample}"
    threads: 1
    log:
        "featureCounts/log/{sample}.log"
    params:
        R2_info
    shell:
        '''
        featureCounts.sh {input} {gdir} {threads} {params} {tmpdir}>&{log}
        '''
rule rRNA:
    input:
        fastq_info
    output:
        "rRNA/{sample}.txt"
    threads: 6
    shell:
        '''
        gdir_root=$(dirname {gdir})
        rRNA=$(species.sh {gdir})_rRNA
        gref=$gdir_root/misc_data/$rRNA/$rRNA
        out_tmp=rRNA/{wildcards.sample}_tmp.txt
        bowtie2.sh -d rRNA $gref {threads} {input} >& $out_tmp
        mv $out_tmp {output}
        '''
rule phix:
    input:
        fastq_info
    output:
        "phix/{sample}.txt"
    threads: 6
    shell:
        '''
        gdir_root=$(dirname {gdir})
        gref=$gdir_root/misc_data/phix/phix
        out_tmp=phix/{wildcards.sample}_tmp.txt
        bowtie2.sh -d phix $gref {threads} {input} >& $out_tmp
        mv $out_tmp {output}
        '''
rule globin:
    input:
        fastq_info
    output:
        "globin/{sample}.txt"
    threads: 6
    shell:
        '''
        gdir_root=$(dirname {gdir})
        globin=$(species.sh {gdir})_globin
        gref=$gdir_root/misc_data/$globin/$globin
        out_tmp=globin/{wildcards.sample}_tmp.txt
        bowtie2.sh -d globin $gref {threads} {input} >& $out_tmp
        mv $out_tmp {output}
        '''
rule qc53:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "qc53/{sample}.RNA_Metrics"
    log:
        "qc53/log/{sample}.log"
    shell:
        '''
        qc53.sh {input} {gdir} 
        '''
rule mark_dup:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "mark_dup/{sample}.dup_metrics"
    log:
        "mark_dup/log/{sample}.log"
    shell:
        '''
        mark_dup.sh {input} >&{log}
        '''
rule featureCounts_all:
    input:
        expand("featureCounts/{sample}",sample=samples),
        "samples"
    output:
        "log/OK.featureCounts"
    log:
        "log/featureCounts.log"
    shell:
        '''
        cat samples | awk '{{print "featureCounts/"$0"\t"$0}}' >.samples_feature
        row_paste.awk infoid=1 colid=0 skip=1 <.samples_feature >featureCounts.txt 2>{log}
        echo "Finished featureCounts" >{output}
        '''
        
rule star_align_all:
    input:
        expand("star_align/{sample}/Log.final.out",sample=samples),
        "samples"
    output:
        "log/OK.star_align"
    log:
        "log/star_align.log"
    shell:
        '''
        cat samples | awk '{{print "star_align/"$0"/Log.final.out\t"$0}}' |
            row_paste.awk infoid=1 colid=0 skip=-1 >star_align/star_QC.txt 2>{log}
        echo "OK" > {output}
        '''
def fastqc_all_input(wildcards):
    files=expand("fastqc/{sample}_fastqc.html",sample=fastqc_samples)
    if fastq_ini=="fastq_raw/":
        files.extend(expand("fastqc_raw/{sample}_fastqc.html",sample=fastqc_samples))
    return(files)    
                
rule fastqc_all:
    input:
        fastqc_all_input
    output:
        "log/OK.fastqc"
    shell:
        '''
        echo OK >{output}
        '''
rule post_align_QC:
    input:
        expand("qc53/{sample}.RNA_Metrics",sample=samples),
        "log/OK.rsem",
        "log/OK.featureCounts",
    output:
        "log/OK.post_align_QC"
    log:
        "log/post_align_QC.log"
    shell:
        '''
        Rscript --vanilla {root}/bin/qc53.R
        mkdir -p multiqc
        multiqc.sh post_align star_align featureCounts rsem >&{log} 
        echo OK >{output}
        '''
        
rule pre_align_QC:
    input:
        "log/OK.fastqc"
    output:
        "log/OK.pre_align_QC"
    log:
        "log/pre_align_QC.log"
    shell:
        '''
        mkdir -p multiqc
        fastqc_raw=""
        if [[ {fastq_ini} == "fastq_raw/" ]]; then
           fastqc_raw="fastqc_raw fastq_trim"
        fi
        multiqc.sh pre_align fastqc $fastqc_raw >&{log}
        echo OK >{output}
        '''
dup_input=expand("mark_dup/{sample}.dup_metrics",sample=samples)
if I==len(samples):
    dup_input.extend(expand("star_align/{sample}/{sample}_dup_log.txt",sample=samples))
    

rule qc_all:
    input:
        expand("globin/{sample}.txt",sample=samples),
        expand("phix/{sample}.txt",sample=samples),
        expand("rRNA/{sample}.txt",sample=samples),
        expand("star_align/{sample}/chr_info.txt",sample=samples),
        dup_input,
        "log/OK.pre_align_QC",
        "log/OK.post_align_QC"
    output:
        "log/OK.qc_all"
    shell:
        '''
        Rscript --vanilla {root}/bin/qc.R
        echo OK>{output}
        '''
        
rule rsem_all:
    input:
        expand("rsem/log/OK.{sample}",sample=samples),
        "samples"
    output:
        "log/OK.rsem"
    log:
        "log/rsem.log"
    shell:
        '''
        cat samples | awk '{{print "rsem/"$0".genes.results\t"$0}}'  > .samples.rsem
        row_paste.awk colid=6 < .samples.rsem >rsem_genes_tpm.txt 2>{log}
        row_paste.awk colid=7 < .samples.rsem >rsem_genes_fpkm.txt 2>>{log}
        row_paste.awk colid=5 < .samples.rsem >rsem_genes_count.txt 2>>{log}
        echo "Finished rsem">{output}
        '''
