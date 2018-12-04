#generate a small data for testing
zcat ../../hg38_gencode_v29/source/gencode.v29.primary_assembly.annotation.gtf.gz |awk '$1=="chr22" || /^#/' |gzip -c >chr22.gft.gz

zcat ../../hg38_gencode_v29/source/GRCh38.primary_assembly.genome.fa.gz |awk '/^>/ { if(go==1 && $1!=">chr22") exit; if(go=$1==">chr22") print; }; !/^>/ && go==1' |gzip -c >chr22.fa.gz
