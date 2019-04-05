#!/bin/bash 
set -eux -o pipefail

outdir=$1
indir=${outdir}/fastq_raw

# export environmental variables
conda=/labs/smontgom/nicolerg/MOTRPAC/PIPELINES/conda
refdata=/labs/smontgom/nicolerg/MOTRPAC/RNA/REFERENCES
export $(/labs/smontgom/nicolerg/MOTRPAC/PIPELINES/motrpac_rnaseq/bin/load_motrpac.sh -c $conda -r $refdata)
export PYTHONPATH=''

# sbatch option		JSON cluster parameter
# --cpus-per-task	"nCPUs"
# --time		"time"
# --mem			"mem"
# --account		"account"
# --mail-type		"mail"
# --output		"output"
# --error		"error"

# move to the folder where we want outputs
cd ${outdir}

# make a symlink to fastq_raw in CWD called fastq_raw
ln -sf ${indir} fastq_raw

jsonfile=$MOTRPAC_root/config/slurm.json

mkdir -p log/cluster

snakemake -j 999 --snakefile $MOTRPAC_root/rna-seq.snakefile \
					--config genome=rn6_ensembl_r95 \
					--cluster-config $jsonfile \
					--cluster \
					"sbatch --account={cluster.account} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--output={cluster.output} \
						--error={cluster.error} \
						--mail-type={cluster.mail}"
															

