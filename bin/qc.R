#!/usr/bin/env Rscript
#--vanilla is not working
##working on collect all QC in one single file
##replace the rRNA in qc53.txt
options(stringsAsFactors=FALSE)
readqcinfo<-function(type, name)
{
    read.delim(paste0("multiqc/",type,"_data/multiqc_",name,".txt"),
               sep="\t",head=TRUE,row.names=1,check.names=FALSE,strip.white=TRUE)
}
checknames<-function(x)
{
    sum(rownames(x)!=samples)!=0
}
#read star qc info the multiqc misses chimeric%, go with manully collected data from Log.final.out 
star<-t(read.delim("star_align/star_QC.txt",sep="\t",row.names=1,strip.white=TRUE,check.names=FALSE))
star<-star[,c(6,7,9:17,25,27,29,30,31,34)]
colnames(star)<-c("reads","avg_input_read_length","uniquely_mapped","%uniquely_mapped","avg_mapped_read_length",
                  "num_splices","num_annotated_splices","num_GTAG_splices","num_GCAG_splice","num_ATAC_splices",
                  "num_noncanonical_splices","%multimapped","%multimapped_toomany","%unmapped_mismatches","%unmapped_tooshort",
                  "%unmapped_other","%chimeric")
star<-sub("%$","",star)
rownames(star)<-sub(" \\|$","",sub("^star_align \\| ","",rownames(star)))
samples<-sort(rownames(star))
star<-star[samples,]
NS<-length(samples)
#read fastqc info
fastqc<-readqcinfo("pre_align","fastqc")[,c("Total Sequences","%GC","total_deduplicated_percentage")]
if(2*length(grep("_R2$",rownames(fastqc)))==nrow(fastqc)){
    ##paired, so we do the average of R1 and R2
    id<-(1:(nrow(fastqc)/2))*2
    if(sum(sub("_R2$","_R1",rownames(fastqc)[id])
           !=rownames(fastqc)[id-1])!=0){
        stop("the R1 and R2 files do not match in the fastqc info")
    }
    fastqc<-(fastqc[id-1,]+fastqc[id,])/2
}
rownames(fastqc)<-sub("_R1$","",rownames(fastqc))

#read the trim info and also clean-up fastqc when necessary
TRIM<-dir.exists("fastq_trim")
if(TRIM){
    trim<-readqcinfo("pre_align","cutadapt")
    rownames(trim)<-sub("_R[12]$","",rownames(trim))
    if(checknames(trim)){
        stop("the trim info is wrong")
    }
    fastqc_raw<-fastqc[grep("^fastqc_raw \\| ",rownames(fastqc)),]
    rownames(fastqc_raw)<-sub("^fastqc_raw \\| ","",rownames(fastqc_raw))
    if(checknames(fastqc_raw)){
        stop("the fastqc_raw info is wrong")
    }
    fastqc<-fastqc[grep("^fastqc \\| ",rownames(fastqc)),]
}
rownames(fastqc)<-sub("^fastqc \\| ","",rownames(fastqc))
if(checknames(fastqc)){
        stop("the fastqc info is wrong")
    }
##read qc53,note that the multiqc is not working well for collectRNAmetrics
qc53<-t(read.delim("qc53.txt",sep="\t",head=TRUE,row.names=1,check.names=FALSE))
qc53<-qc53[,-match(c("RIBOSOMAL_BASES","PCT_RIBOSOMAL_BASES","SAMPLE","LIBRARY","READ_GROUP"),
                  colnames(qc53)),drop=FALSE]
Nqc53<-ncol(qc53)
nqc53<-nrow(qc53)
if(nqc53==NS){
    UMI=FALSE
    qc53<-qc53[samples,,drop=FALSE]
}else if(nqc53==2*NS){
    UMI=TRUE
    qc53u<-qc53[paste0("UMI_",samples),,drop=FALSE]
    colnames(qc53u)<-paste0(colnames(qc53),"_umi")
    rownames(qc53u)<-sub("^UMI_","",rownames(qc53u))
    qc53<-cbind(qc53[samples,,drop=FALSE],qc53u)
}else{
    stop("star_qc and qc53 do not match")
}

id<-c("CODING","UTR","INTRONIC","INTERGENIC","MRNA")
loc<-match(c(paste0("PCT_",id,"_BASES"),"MEDIAN_5PRIME_TO_3PRIME_BIAS"),
           colnames(qc53)[1:Nqc53])
if (UMI){
    loc<-c(loc,loc+Nqc53)
}
qc53<-qc53[,loc,drop=FALSE]
id2<-tolower(id)
colnames(qc53)[1:6]<-c(paste0("%",id2),"median_5'_3'_bias")
if(UMI){
    colnames(qc53)[1:6+6]<-
        c(paste0("%","umi_",id2),"umi_median_5'_3'_bias")
}
id<-grep("^%",colnames(qc53))
qc53[,id]<-round(qc53[,id]*100,dig=2)

##collect data for rRNA, globin and phix and ERCC and duplicates when present
misc<-matrix(NA,NS,6)
colnames(misc)<-c("globin","rRNA","phix","ERCC","picard_dup","UMI_dup")
for(i in 1:NS){
    SID<-samples[i]
    for(j in 1:4){ #change to 4 when good
        zz<-pipe(paste0("tail -1 ",colnames(misc)[j],"/",SID,".txt |tr -d %"))
        misc[i,j]<-scan(zz,n=1,quiet=TRUE)
        close(zz)
    }
    misc[i,5]<-read.delim(paste0("mark_dup/",SID,".dup_metrics"),skip=6,head=TRUE,row=1)[1,"PERCENT_DUPLICATION"]*100
    if(UMI){
        zz<-pipe(paste0("tail -1 star_align/UMI_",SID,"/",SID,"_dup_log.txt |cut -f 6"))
        misc[i,6]<-scan(zz,n=1,quiet=TRUE)*100
        close(zz)
    }
}
misc<-round(misc,dig=2)
colnames(misc)<-paste0("%",colnames(misc))

##Read the chr_info.txt
chr_info<-matrix(NA,NS,10)
colnames(chr_info)<-c("chrX","chrY","chrM","chrAuto","contig",
                      "umi_chrX","umi_chrY","umi_chrM","umi_chrAuto","umi_contig")
readchr<-function(sid){
    x<-read.delim(paste0("star_align/",sid,"/chr_info.txt"),
                  sep="\t",row.names=1,header=FALSE)
    xt<-sum(x[,2])
    y<-x[c("chrX","chrY","chrM"),2]/xt*100
    y<-c(y,sum(x[grep("^chr[1-9]",rownames(x)),2])/xt*100)
    c(y,100-sum(y))
}
    
for(i in 1:NS){
    chr_info[i,1:5]<-readchr(samples[i])
    if(UMI){
        chr_info[i,6:10]<-readchr(paste0("UMI_",samples[i]))
    }
}
colnames(chr_info)<-paste0("%",colnames(chr_info))
if(!UMI){
    chr_info<-chr_info[,1:5] #the other columnes are not needed
}
chr_info<-round(chr_info,dig=2)
#Now putting all togther
qc<-NULL
if(TRIM) qc<-cbind("reads_raw"=fastqc[,1],"%trimmed"=round(as.numeric(star[,1])/fastqc[,"Total Sequences"]*100,dig=2),
                   "%trimmed_bases"=trim[,"percent_trimmed"])
qc<-cbind(qc,reads=star[,1],"%GC"=round(fastqc[,"%GC"],dig=2),"%dup_sequence"=100-round(fastqc[,"total_deduplicated_percentage"],dig=2),
          misc[,1:5])
if(UMI){
    qc<-cbind(qc,"%umi_dup"=misc[,6])
}
qc<-cbind(qc,star[,-1],chr_info,qc53)

id<-grep("_percent$",colnames(qc))
colnames(qc)[id]<-paste0("%",sub("_percent$","",colnames(qc)[id]))

write.table(qc,"qc_info.txt",row.names=TRUE,col.names=NA,quote=FALSE,sep="\t")
