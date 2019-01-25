##This file needs to be run interactively to spot problems
##Download the files specified in the human and rat
#as of Nov 20, 2018

#The following is for the hg38 gencode version 29

#MOTRPAC_root should be set to the root of the code base
#MOTRPAC_refdata should be set to the intended genome reference data

#The environments should be setup accorrding to load_motrpac.sh and after the softwares have been installed by conda_install.sh

cd $MOTRPAC_refdata
mkdir -p hg38_gencode_v29/source
cd hg38_gencode_v29/source
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz

#The following is for the rn6 ensembl release version 95
cd $MOTRPAC_refdata
mkdir -p rn6_ensembl_r95/source
cd rn6_ensembl_r95/source
wget ftp://ftp.ensembl.org/pub/release-95/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-95/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.95.gtf.gz

#The following is for the mm10 gencode version v19, not used for MotrPAC, but for other purpose
cd $MOTRPAC_refdata
mkdir -p mm10_gencode_v19/source
cd mm10_gencode_v19/source
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/GRCm38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.primary_assembly.annotation.gtf.gz

#copy rRNA and globin data from code base into refdata
mkdir -p $MOTRPAC_refdata/misc_data/log
#each fa file is one folder
cd $MOTRPAC_refdata/misc_data
for fa in $MOTRPAC_root/misc_data/*.fa;
do
    faid=$(basename $fa .fa)
    mkdir -p $faid
    cp $fa $faid
    cd $faid
    bowtie2-build $faid.fa $faid >& ../log/${faid}_bowtie_index.log
    bismark_genome_preparation . >& ../log/${faid}_bismark_index.log
    cd ..
done
