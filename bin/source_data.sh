##This file needs to be run interactively to spot problems
##Download the files specified in the human and rat
#as of Oct 18, 2021. Updated the human from v30 to v38

#The following is for the hg38 gencode version 38

#MOTRPAC_root should be set to the root of the code base
#MOTRPAC_refdata should be set to the intended genome reference data

#The environments should be setup accorrding to load_motrpac.sh and after the softwares have been installed by conda_install.sh

cd $MOTRPAC_refdata
mkdir -p hg38_gencode_v38/source
cd hg38_gencode_v38/source
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz

#The following is for the rn6 ensembl release version 96
cd $MOTRPAC_refdata
mkdir -p rn6_ensembl_r96/source
cd rn6_ensembl_r96/source
wget ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-96/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.96.gtf.gz

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
