# motrpac
This contains the pipelines that are useful for the motrpac projects

## The conda enviroment
- installatin folder: /sc/orga/projects/sealfs01a/conda
- The bin path for snakemake and others: /sc/orga/projects/sealfs01a/conda/python3/bin

## How to use
- make sure that snakemake is in your path `export PATH=/sc/orga/projects/sealfs01a/conda/python3/bin:$PATH`
- go to the folder that contains a subfolde `fastq` where the fastq file names are specified according to the pipeline
- snakemake -f path/rna-seq.snakefile
