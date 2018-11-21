localrules: all
rule all:
    input:
        "log/OK.star_index",
        "log/OK.rsem_index",
        "log/OK.bowtie_index",
        "log/OK.qc53_ref"
        
samples, = glob_wildcards("{sample}")        
wildcard_constraints:
    sample="[^/]+" #do not allow to go to the sub directory
    
rule genome:
    output:
        "genome.gtf",
        "genome.fa"
    log:
        "log/genome.log"
    shell:
        '''
        genome.sh {log}
        '''
        
rule star_index:
    input:
        "genome.gtf",
        "genome.fa"
    output:
        "log/OK.star_index"
    log:
        "star_index/log/build.log"
    threads: 6
    shell:
        '''
        star_index.sh {threads} {log}
        echo OK >{output}
        '''
        
rule rsem_index:
    input:
        "genome.gtf",
        "genome.fa"
    output:
        "log/OK.rsem_index"
    log:
        "rsem_index/log/build.log"
    #threads: 6,not supported as of rsem 1.3.1
    shell:
        '''
        rsem_index.sh {log}
        echo OK >{output}
        '''
        
rule bowtie_index:
    input:
        "genome.fa"
    output:
        "log/OK.bowtie_index"
    log:
        "bowtie2_index/log/build.log"
    threads: 6
    shell:
        '''
        bowtie2_index.sh {threads} {log}
        echo OK >{output}
        '''
        
rule qc53_ref:
    input:
        "genome.gtf"
    output:
        "log/OK.qc53_ref"
    log:
        "qc53_ref/log/build.log"
    shell:
        '''
        qc_53.ref.sh {log}
        echo OK >{output}
        '''
