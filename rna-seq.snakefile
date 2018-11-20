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

localrules: all,fastqc_all,star_align_all,featureCounts_all,rsem_all,rRNA_all,qc53_all

rule all:
    input:
        "log/OK.fastqc",
        "log/OK.star_align",
        "log/OK.featureCounts",
        "log/OK.rsem",
        "log/OK.rRNA",
        "log/OK.qc53"
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
        qc_log="star_align/{sample}/Log.final_filled.out" 
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
rule qc53:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "qc53/{sample}{grp,_[0-9]+-[0-9]+k}.RNA_Metrics" #grp . _1-2k., etc
    log:
        "qc53/log/{sample}{grp}.log"
    shell:
        '''
        {bin}/qc53.sh {wildcards.sample} {wildcards.grp} {gdir} 
        '''
rule qc53_sample_all:
    input:
        expand("qc53/{{sample}}{grp}.RNA_Metrics",grp=["_0-0k","_0-1k","_1-2k","_2-3k","_3-4k","_4-6k","_6-10k","_10-0k"])
    output:
        "qc53/OK.{sample}"
    shell:
        '''
        echo "finished sample {wildcards.sample}" > {output}
        ''' 
        
rule qc53_all:
    input:
        expand("qc53/{sample}{grp}.RNA_Metrics",sample=samples,grp=["_0-0k","_0-1k","_1-2k","_2-3k","_3-4k","_4-6k","_6-10k","_10-0k"])
    output:
        "log/OK.qc53"
    log:
        "qc53/log/all.log"
    shell:
        '''
        module load R
        echo "finished all qc53" > {output}
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
    log:
        "log/featureCounts.log"
    shell:
        '''
        SID.sh 1 1 <sample_info.txt| awk '{{print "featureCounts/"$0"\t"$0}}' |
        row_paste.awk infoid=1 colid=0 skip=1 >featureCounts.txt 2>{log}
        echo "Finished featureCounts" >{output}
        '''
        
rule star_align_all:
    input:
        expand("star_align/{sample}/Log.final_filled.out",sample=samples)
    output:
        "log/OK.star_align"
    log:
        "log/star_align.log"
    shell:
        '''
        SID.sh 1 1 <sample_info.txt| awk '{{print "star_align/"$0"/Log.final_filled.out\t"$0}}' |
        row_paste.awk infoid=1 colid=0 skip=1 >star_align/star_QC.txt 2>{log}
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
        SID.sh 1 1 <sample_info.txt | awk '{{print "rsem/"$0".genes.results\t"$0}}'  > .samples.rsem
        row_paste.awk colid=6 < .samples.rsem >rsem_genes_tpm.txt 2>{log}
        row_paste.awk colid=7 < .samples.rsem >rsem_genes_fpkm.txt 2>>{log}
        row_paste.awk colid=5 < .samples.rsem >{output} 2>>{log}
        echo "Finished rsem">{output}
        '''
        
