#Author: Yongchao Ge
#July 20,2018: The initial push to sealfonlab
#The usage are in the file rna-seq_README.md

##Things to do:
#1 cutadapt is missing for the moment as it requires adapter file and it required more tweaking
#2 fastqc currently doesn't make use of the adpater
#3 multiqc not implemented yet

configfile: "config.yaml"
bin=os.environ['MOTRPAC_ROOT']+"/bin"
def read_samples(file):
    samples={}
    fastqc_samples=[]
    n=0
    with open(file,'rt') as f:
        for line in f:
            n+=1
            if n==1:
                continue
            info=line.split()
            sid=info[0]
            R=info[1]
            if R!="R1" and R!="R1,R2":
                print("Error: in R column for sample "+sid+" with R="+R)
                exit(1)
            samples[sid]=R.split(",")
            fastqc_samples.extend(["%s_%s" % (sid,x) for x in samples[sid]])         
    return samples,fastqc_samples

def fastq_info(wildcards):
    sid=wildcards.sample
    return(expand("fastq/{sample}_{R}.fastq.gz",sample=sid,R=samples[sid]))
           
gdir=config["genomedir"]
samples,fastqc_samples=read_samples("sample_info.txt")

localrules: all,fastqc_all,star_align_all,featureCounts_all,rsem_all,rRNA_all

rule all:
    input:
        "log/OK.fastqc",
        "log/OK.star_align",
        "log/OK.featureCounts",
        "log/OK.rsem",
        "log/OK.rRNA"
    output:
        "log/OK.all"
    shell:
        '''
        echo "Finished all" >{output}
        '''
rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastqc/{sample}_fastqc.html"
    log:
        "fastqc/log/{sample}.log"
    shell:
        '''
        module load fastqc/0.11.7
        fastqc -o fastqc {input} >& {log}
        '''

rule star_align:
    input:
        fastq_info 
    output:
        bam="star_align/{sample}/Aligned.sortedByCoord.out.bam",
        rsem_bam="star_align/{sample}/Aligned.toTranscriptome.out.bam",
        qc_log="star_align/{sample}/Log.final.out" 
    log:
        "star_align/log/{sample}.log"
    threads: 6
    shell:
        '''
        {bin}/star_align.sh {wildcards.sample} {gdir} {threads} {input} >&{log}
        '''
        
rule rsem:
    input:
        "star_align/{sample}/Aligned.toTranscriptome.out.bam"
    output:
        "rsem/log/OK.{sample}"
    log:
        "rsem/log/{sample}.log"
    threads: 6
    shell:
        '''
        {bin}/rsem.sh {wildcards.sample} {gdir} {threads} >&{log}
        echo "Finished rsem" > {output}
        '''
        
rule featureCounts:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "featureCounts/{sample}"
    log:
        "featureCounts/log/{sample}.log"
    shell:
        '''
        module load subread/1.5.0-p1
        pairopt=""
        if [ -e fastq/{wildcards.sample}_R2.fastq.gz ]; then
            pairopt="-p"
        fi
        featureCounts -a {gdir}/genome.gtf -o {output} $pairopt -M --fraction {input} >& {log}
        '''

rule rRNA:
    input:
        fastq_info
    output:
        "rRNA/{sample}.txt"
    threads: 6
    shell:
        '''
        {bin}/rRNA.sh {wildcards.sample} {gdir} {threads}
        '''
        
rule fastqc_all:
    input:
        expand("fastqc/{sample}_fastqc.html",sample=fastqc_samples)
    output:
         "log/OK.fastqc"
    shell:
        '''
        echo "Finished all fastqc">{output}
        '''
        
rule featureCounts_all:
    input:
        expand("featureCounts/{sample}",sample=samples)
    output:
        "log/OK.featureCounts"
    shell:
        '''
        echo "Finished featureCounts" >{output}
        '''
        
rule star_align_all:
    input:
        expand("star_align/{sample}/Aligned.sortedByCoord.out.bam",sample=samples)
    output:
        "log/OK.star_align"
    shell:
        '''
        echo "Finished star_align" > {output}
        '''
rule rRNA_all:
    input:
        expand("rRNA/{sample}.txt",sample=samples)
    output:
        "log/OK.rRNA"
    shell:
        '''
        echo "Finished all rRNA" >{output}
        '''

rule rsem_all:
    input:
        expand("rsem/log/OK.{sample}",sample=samples)
    output:
        "log/OK.rsem"
    log:
        "log/rsem.log"
    shell:
        '''
        echo "Finished rsem">{output}
        '''
        
