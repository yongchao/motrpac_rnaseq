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
checknames<-function(x,title)
{
    if(sum(rownames(x)!=samples)==0) return(1)
    #to fix the problem of multqc that the order of S1 is later than S10
    if(sum(sort(rownames(x))!=sort(samples))==0) return(2)
    stop(paste0(title," is wrong"))
}
    
#read star qc info the multiqc misses chimeric%, go with manully collected data from Log.final.out 
star<-t(read.delim("star_align/star_QC.txt",sep="\t",row.names=1,strip.white=TRUE,check.names=FALSE))
star<-star[,c(6,7,9:17,25,27,29,30,31,34)]
colnames(star)<-c("reads","avg_input_read_length","uniquely_mapped","%uniquely_mapped","avg_mapped_read_length",
                  "num_splices","num_annotated_splices","num_GTAG_splices","num_GCAG_splices","num_ATAC_splices",
                  "num_noncanonical_splices","%multimapped","%multimapped_toomany","%unmapped_mismatches","%unmapped_tooshort",
                  "%unmapped_other","%chimeric")
star<-sub("%$","",star)
rownames(star)<-sub(" \\|$","",sub("^star_align \\| ","",rownames(star)))

samples<-sort(rownames(star))

star<-star[samples,]
NS<-length(samples)
#read fastqc info
fastqc<-readqcinfo("pre_align","fastqc")[,c("Total Sequences","%GC","total_deduplicated_percentage")]
All_Single<-TRUE
if(length(grep("_R2$",rownames(fastqc)))>0){
    #some are paired
    All_single<-FALSE
    ##paired, so we do the average of R1 and R2
    #id<-(1:(nrow(fastqc)/2))*2
    #if(sum(sub("_R2$","_R1",rownames(fastqc)[id])
    #       !=rownames(fastqc)[id-1])!=0){
    #    stop("the R1 and R2 files do not match in the fastqc info")
    #}
    ##fastqc<-(fastqc[id-1,]+fastqc[id,])/2
    lab<-sub("_R2$","_R1",rownames(fastqc))
    fastqc<-apply(fastqc,2,function(x){tapply(x,lab,mean)})
}
rownames(fastqc)<-sub("_R1$","",rownames(fastqc))

#read the trim info and also clean-up fastqc when necessary
TRIM<-dir.exists("fastq_trim")
UMI<-dir.exists("fastq_attach")
if(TRIM){
    trim<-readqcinfo("pre_align","cutadapt")
    rownames(trim)<-sub("_R[12]$","",rownames(trim))
    if(checknames(trim,"trim")==2){
        trim<-trim[samples,]
    }
    fastqc_raw<-fastqc[grep("^fastqc_raw \\| ",rownames(fastqc)),]
    rownames(fastqc_raw)<-sub("^fastqc_raw \\| ","",rownames(fastqc_raw))
    if(checknames(fastqc_raw,"fastq_raw")==2){
        fastqc_raw<-fastqc_raw[samples,]
    }
    fastqc<-fastqc[grep("^fastqc \\| ",rownames(fastqc)),]
}
rownames(fastqc)<-sub("^fastqc \\| ","",rownames(fastqc))
if(checknames(fastqc,"fastqc")==2){
    fastqc<-fastqc[samples,]
}

##read qc53,note that the multiqc is not working well for collectRNAmetrics
qc53<-t(read.delim("qc53.txt",sep="\t",head=TRUE,row.names=1,check.names=FALSE))
qc53<-qc53[,-match(c("RIBOSOMAL_BASES","PCT_RIBOSOMAL_BASES","SAMPLE","LIBRARY","READ_GROUP"),
                  colnames(qc53)),drop=FALSE]
Nqc53<-ncol(qc53)
nqc53<-nrow(qc53)
if(nqc53==NS){
    qc53<-qc53[samples,,drop=FALSE]
}else{
    stop("star_qc and qc53 do not match")
}

id<-c("CODING","UTR","INTRONIC","INTERGENIC","MRNA")
loc<-match(c(paste0("PCT_",id,"_BASES"),"MEDIAN_5PRIME_TO_3PRIME_BIAS"),
           colnames(qc53)[1:Nqc53])
qc53<-qc53[,loc,drop=FALSE]
qc53<-round(qc53,dig=2)
id2<-tolower(id)
colnames(qc53)[1:6]<-c(paste0("%",id2),"median_5_3_bias")

id<-grep("^%",colnames(qc53))
qc53[,id]<-round(qc53[,id]*100,dig=2)

##collect data for rRNA, globin and phix and duplicates when present
misc<-matrix(NA,NS,6)
colnames(misc)<-c("globin","rRNA","phix","picard_dup","UMI_dup","adapter_detected")
for(i in 1:NS){
    SID<-samples[i]
    for(j in 1:3){ #change to 4 when good
        zz<-pipe(paste0("tail -1 ",colnames(misc)[j],"/",SID,".txt |tr -d %"))
        misc[i,j]<-scan(zz,n=1,quiet=TRUE)
        close(zz)
    }
    misc[i,4]<-read.delim(paste0("mark_dup/",SID,".dup_metrics"),skip=6,head=TRUE,row=1)[1,"PERCENT_DUPLICATION"]*100
    if(UMI){
        zz<-pipe(paste0("tail -1 star_align/",SID,"/",SID,"_dup_log.txt |cut -f 6"))
        misc[i,5]<-scan(zz,n=1,quiet=TRUE)*100
        close(zz)
    }
    if(TRIM){
        zz<-pipe(paste0("grep \"with adapter\" fastq_trim/log/log.",SID,"|awk -F '[(%]' '{print $2}'"))
        zval<-scan(zz,quiet=TRUE)
        ##future plan is to consider the paired info for each indivudal files
        #if(length(zval)!= PAIRED+1){
        #    stop("the fastq_trim log for the contained adapter% is not consistent with the pairedness of the fastq data for sample", SID)
        #}
        misc[i,6]<-mean(zval)
        close(zz)
    }
}
misc<-round(misc,dig=2)
colnames(misc)<-paste0("%",colnames(misc))

##Read the chr_info.txt
chr_info<-matrix(NA,NS,5)
colnames(chr_info)<-c("chrX","chrY","chrM","chrAuto","contig")

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
}
colnames(chr_info)<-paste0("%",colnames(chr_info))
chr_info<-round(chr_info,dig=2)
#Now putting all togther
qc<-NULL
if(TRIM) qc<-cbind("reads_raw"=fastqc_raw[,1],"%adapter_detected"=misc[,6],"%trimmed"=round(100-as.numeric(star[,1])/fastqc_raw[,"Total Sequences"]*100,dig=2),
                   "%trimmed_bases"=round(trim[,"percent_trimmed"],dig=2))
qc<-cbind(qc,reads=star[,1],"%GC"=round(fastqc[,"%GC"],dig=2),"%dup_sequence"=round(100-fastqc[,"total_deduplicated_percentage"],dig=2),
          misc[,1:4])
if(UMI){
    qc<-cbind(qc,"%umi_dup"=misc[,5])
}
qc<-cbind(qc,star[,-1],chr_info,qc53)

id<-grep("_percent$",colnames(qc))
colnames(qc)[id]<-paste0("%",sub("_percent$","",colnames(qc)[id]))
colnames(qc)<-sub("^%","pct_",colnames(qc))

#We could have write a date and the project folder here
write.table(qc,"qc_info.csv",row.names=TRUE,col.names=NA,quote=FALSE,sep=",")
