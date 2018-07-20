# motrpac
This repository contains the pipelines that are useful for the motrpac projects

##version
-  See the file [version](version), please update it when we push the changes to the git.
## code folder structure
- `bin`: awk, bash, R or python scripts that are used in the pipeline
- `config`: different configuration files
- `refdata`: different genome and gtf and other indexes (star,rsem, etc). This only sits on the frozen project folder
- snakefiles and the corresponding readme files in the main folder: The implementation of different pipelines

## The conda enviroment
- installation folder: /sc/orga/projects/sealfs01a/conda
- The bin path for snakemake and others: `/sc/orga/projects/sealfs01a/conda/python3/bin`

## The frozen database and github code
- A frozen copy of github code will be living at `/sc/orga/projects/sealfs01a/motrpac/`
- The reference databases that are using the RNA-seq or RBBS or other MOP (including the genome and gtf version and software version) will be living at `/sc/orga/projects/sealfs01a/motrpac/refdata`

## All users can modify the code in the github. For the moment, the conda and motrpac folders are only read only for the group. Yongchao is now taking the responsibility of 
changing the the conda enviroment and the frozen database and github code. This arrangement can be changed later.

## How to use
- make sure that snakemake and other required scripts are in your path and setting up the envivonment variables 'MOTRPAC_ROOT' and 'CONDA_ROOT' by running path/bin/load_motrpac.sh
- copy the commands emitted from the previous command to the terminal (this can be done better for future)
- go to a project folder that contains a subfolder `fastq` where the fastq file names are specified according to the pipeline
- Run the command according to the specific pipeline. For example, the rna-seq usage can be seen in the file [rna-seq_README.md](rna-seq_README.md)
- When necessary, the motrpac environment can be unloaded by running the command path/bin/unload_motrpac.sh and copy the commands emitted to the terminal
