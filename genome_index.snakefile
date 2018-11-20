localrules: all
rule all:
    input:
        "log/OK.star_index",
        "log/OK.rsem_index",
        "log/OK.bowtie_index",
        "log/OK.qc53_ref",
        "log/OK.kallisto_index"
    output:
        "log/OK.all"
    shell:
        '''
        echo OK >{output}
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
    params:
        " --runMode genomeGenerate"
        " --genomeDir star_index"
        " --genomeFastaFiles genome.fa"
        " --sjdbGTFfile genome.gtf"
        " --sjdbOverhang 99"
        " --outFileNamePrefix star_index/"
    shell:
        '''
        module load star/2.5.4b
        STAR --runThreadN {threads} {params} >&{log}
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
    #threads: 6
    shell:
        '''
        #module load rsem/1.2.21, not in the system and using the version installed at conda
        #threads do not seem work for the index
        rsem-prepare-reference --gtf genome.gtf genome.fa rsem_index/genome >&{log}
        echo OK >{output}
        '''
rule bowtie_index:
    input:
        "genome.fa"
    output:
        "log/OK.bowtie_index"
    log:
        "bowtie2_index/log/build.log"
    #threads: 1
        #the multi-core is not working,see https://github.com/BenLangmead/bowtie2/issues/36
        #unless we move to bowtie 2.2.9
    shell:
        '''
        #module load bowtie2/2.2.8, not in the same system and using the version installed at conda
        bowtie2-build genome.fa bowtie2_index/genome >&{log}
        echo OK >{output}
        '''
rule kallisto_index:
    input:
        "genome.fa"
    output:
        "log/OK.kallisto_index"
    log:
        "kallisto_index/log/build.log"
    shell:
        '''
        module load kallisto/0.43.0
        kallisto index -i kallisto_index/genome genome.fa >&{log}
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
        module load ucsc-utils/2015-04-07
        gtfToGenePred -ignoreGroupsWithoutExons -genePredExt genome.gtf genome_pred.txt
        awk 'BEGIN{{FS=OFS="\\t"}};{{print $12, $1, $2,$3,$4,$5,$6,$7,$8,$9,$10}}' genome_pred.txt>genome_flat.txt
        bb=(0 0 1 2 3 4 6 10)
        be=(0 1 2 3 4 6 10 0)
        for i in "${{!bb[@]}}";
        do
	   exon_length.awk lmin=$((1000*bb[i])) lmax=$((1000*be[i])) genome_flat.txt >qc53_ref/genome_flat_${{bb[i]}}-${{be[i]}}k.txt
        done
        echo OK >{output}
        '''
