#Author: Yongchao Ge

#Nov 20, 2018: The intial version that is shared with Stanford
#July 20,2018: The initial push to sealfonlab

#The usage are in the file rna-seq_README.md

#Folder structures at root of working folder
#fastq_raw: raw fastq files, with no adpaters removed, probably a softlink to the fastq files in the output folder of bcl2fastq

#For each sample, we have three possible fastq files
#${sid}_R1.fastq.gz, required for all data. 
#${sid}_R2.fastq.gz, required for pairied end reads
#${sid}_I1.fastq.gz, required for NuGEN with UMI for UMI processsing

#configure genome by using --config genome=hg38_gencode_v29 etc
#Updating the link for hg or hg38 so that we don't have to specify the whole name

root=os.environ['MOTRPAC_root']
if "genome" in config:
    gdir=os.environ['MOTRPAC_refdata']+"/"+config["genome"]
else:
    gdir=os.environ['MOTRPAC_refdata']+"/hg38_gencode_v29"

#If softlins fastq_raw to fastq folder, then no trim is happening
    
#do not allow to go to the sub directory
wildcard_constraints:
    sample="[^/]+"
    
#If fastq_raw exists, then trim and all folders are based on fastq_trim
#Otherwise lookinf for files in fastq folder
import os
if os.path.isdir("fastq"):
    fastq_ini="fastq/"
    fastq="fastq/"
elif os.path.isdir("fastq_raw"):
    fastq_ini="fastq_raw/"
    fastq="fastq_trim/"
else:
    print("The sub folders fastq_raw or fastq do not exit, exit\n")
    sys.exit(1)


samples,=glob_wildcards(fastq_ini+"{sample,[^/]+}_R1.fastq.gz")
fastqc_samples,=glob_wildcards(fastq_ini+"{sample,[^/]+_R[12]}.fastq.gz")

#construct R2 and I1 info
R2={}
R={}
I=0    
for s in samples:
    if os.path.isfile(fastq_ini+s+"_R2.fastq.gz"):
        R2[s]=1
        R[s]=["R1","R2"]
    else:
        R2[s]=0
        R[s]=["R1"]
        
    if os.path.isfile(fastq_ini+s+"_I1.fastq.gz"):
        I+=1

samples_all=[s for s in samples] #this is for summarization, it may include UMI
if I>0:
    if I!=len(samples):
        print("Not all samples have UMI index files, exit\n")
        print(samples)
        sys.exit(1)
    else:
        trim_input="fastq_attach/"
        for s in samples:
            samples_all.extend(["UMI_"+s])
            R2["UMI_"+s]=R2[s] #not necessary, just to avoid keyerror
            R["UMI_"+s]=R[s]
        
else:
    trim_input="fastq_raw/"

    
#All samples should have I or not
#We could mix smaples with single ends and paired ends
        
if(len(samples)==0):
    print("There are no fastq files in fastq_raw or fastq folders, exit\n")
    sys.exit(1)
        
def R2_info(wildcards):
    return(R2[wildcards.sample])
    
def fastq_info(wildcards,prefix=fastq):
    sid=wildcards.sample
    return(expand(prefix+"{sample}_{R}.fastq.gz",sample=sid,R=R[sid]))

ruleorder: trim > trim_single
#ruleorder: UMI > star_align
#ruleorder: UMI > rsem
#ruleorder: UMI > featureCounts
#ruleorder: UMI > qc53

localrules: all,star_align_all,featureCounts_all,rsem_all,samples_all,fastqc_all

rule all:
    input:
        "log/OK.pre_align_QC",
        "log/OK.post_align_QC",
        "log/OK.qc_all",
        "log/OK.star_align",
        "log/OK.featureCounts",
        "log/OK.rsem"

rule samples_all:
    output:
        "samples_all",
        "samples"
    run:
        with open(output[0], "w") as out:
            for s in samples_all:
                out.write(s+"\n")
        out.close()        
        with open(output[1], "w") as out:
            for s in samples:
                out.write(s+"\n")
        out.close()
                
rule fastqc: #on the final fastq data before star
    input:
        fastq+"{sample}.fastq.gz"
    output:
        "fastqc/{sample}_fastqc.html"
    log:
        "fastqc/log/{sample}.log"
    shell:
        '''
        fastqc.sh {input} fastqc >& {log}
        '''
        
rule fastqc_raw: #on the raw data
    input:
        "fastq_raw/{sample}.fastq.gz"
    output:
        "fastqc_raw/{sample}_fastqc.html"
    log:
        "fastqc_raw/log/{sample}.log"
    shell:
        '''
        fastqc.sh {input} fastqc_raw >& {log}
        '''
        
rule UMI_attach:
    input:
        "fastq_raw/{sample}.fastq.gz"
    output:
        "fastq_attach/{sample}.fastq.gz"
    log:
        "fastq_attach/log/{sample}.log"
    shell:
        '''
        SID={wildcards.sample}
        I=${{SID%R[12]}}I1  #remove R[12] and then add back to I1
        zcat {input} | UMI_attach.awk -v Ifq=fastq_raw/$I.fastq.gz|gzip -c>{output} 
        '''
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
        star_align.sh {gdir} {threads} {input} >&{log}
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
def UMI_input(wildcards):
    SID=wildcards.sample
    SID=SID[4:]
    return(expand("star_align/"+SID+"/Aligned.{type}.out.bam",type=["sortedByCoord","toTranscriptome"]))
    
rule UMI:
    input:UMI_input #The naivway is not going to work by appending UMI
    output:"star_align/{sample}/Aligned.sortedByCoord.out.bam",
           rsem_bam=temp("star_align/{sample}/Aligned.toTranscriptome.out.bam")
    log:"star_align/{sample}/dedup.log"
    wildcard_constraints:
        sample="UMI_[^/]+"
    params:
        R2_info
    shell:
        '''
        UMI.sh {input} {params} >&{log}
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
        rsem.sh {input} {gdir} {threads} {params} >&{log}
        echo "Finished rsem" > {output}
        '''
        
rule featureCounts:
    input:
        "star_align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "featureCounts/{sample}"
    log:
        "featureCounts/log/{sample}.log"
    params:
        R2_info
    shell:
        '''
        featureCounts.sh {input} {gdir} {threads} {params} >&{log}
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
        bowtie2.sh $gref {threads} {input} >& $out_tmp
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
        bowtie2.sh $gref {threads} {input} >& $out_tmp
        mv $out_tmp {output}
        '''
rule ERCC:
    input:
        fastq_info
    output:
        "ERCC/{sample}.txt"
    threads: 6
    shell:
        '''
        gdir_root=$(dirname {gdir})
        gref=$gdir_root/misc_data/ERCC92/ERCC92
        out_tmp=ERCC/{wildcards.sample}_tmp.txt
        bowtie2.sh $gref {threads} {input} >& $out_tmp
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
        bowtie2.sh $gref {threads} {input} >& $out_tmp
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
        expand("featureCounts/{sample}",sample=samples_all),
        "samples_all"
    output:
        "log/OK.featureCounts"
    log:
        "log/featureCounts.log"
    shell:
        '''
        cat samples_all| awk '{{print "featureCounts/"$0"\t"$0}}' |
        row_paste.awk infoid=1 colid=0 skip=1 >featureCounts.txt 2>{log}
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
        expand("qc53/{sample}.RNA_Metrics",sample=samples_all),
        "log/OK.rsem",
        "log/OK.featureCounts",
    output:
        "log/OK.post_align_QC"
    shell:
        '''
        Rscript --vanilla {root}/bin/qc53.R
        mkdir -p multiqc
        multiqc.sh post_align star_align featureCounts rsem 
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
        multiqc.sh pre_align fastqc $fastqc_raw
        echo OK >{output}
        '''
rule qc_all:
    input:
        expand("globin/{sample}.txt",sample=samples),
        expand("ERCC/{sample}.txt",sample=samples),
        expand("phix/{sample}.txt",sample=samples),
        expand("rRNA/{sample}.txt",sample=samples),
        expand("star_align/{sample}/chr_info.txt",sample=samples_all),
        expand("mark_dup/{sample}.dup_metrics",sample=samples),
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
        expand("rsem/log/OK.{sample}",sample=samples_all),
        "samples_all"
    output:
        "log/OK.rsem"
    log:
        "log/rsem.log"
    shell:
        '''
        cat samples_all | awk '{{print "rsem/"$0".genes.results\t"$0}}'  > .samples.rsem
        row_paste.awk colid=6 < .samples.rsem >rsem_genes_tpm.txt 2>{log}
        row_paste.awk colid=7 < .samples.rsem >rsem_genes_fpkm.txt 2>>{log}
        row_paste.awk colid=5 < .samples.rsem >rsem_genes_count.txt 2>>{log}
        echo "Finished rsem">{output}
        '''
