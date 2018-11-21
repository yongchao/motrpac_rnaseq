##This file needs to be run interactively to spot problems
##Download the files specified in the human and rat
#as of Nov 20, 2018

#The following is for the hg38 gencode version 29
cd $MOTRPAC_ROOT/refdata
mkdir -p hg38_gencode_v29/source
cd hg38_gencode_v29/source
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz

#The following is for the rn6 release version 94
cd $MOTRPAC_ROOT/refdata
mkdir -p rn6_ensembl_r94/source
cd rn6_ensembl_r94/source
wget ftp://ftp.ensembl.org/pub/release-94/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-94/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
