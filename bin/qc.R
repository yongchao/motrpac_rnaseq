#!/usr/bin/env Rscript
#--vanilla is not working
##working on collect all QC in one single file
##replace the rRNA in qc53.txt
options(stringsAsFactors=FALSE)
star<-t(read.delim("star_align/star_QC.txt",sep="\t",head=TRUE,row.names=1,check.names=FALSE))
qc53<-t(read.delim("qc53/qc53.txt",sep="\t",head=TRUE,row.names=1,check.names=FALSE))
qc53<-qc53[,-match(c("RIBOSOMAL_BASES","PCT_RIBOSOMAL_BASES","SAMPLE","LIBRARY","READ_GROUP"),
                  colnames(qc53)),drop=FALSE]
samples<-rownames(star)
NS<-length(samples)
Nqc53<-ncol(qc53)
nqc53<-nrow(qc53)
if(nqc53==NS){
    UMI=FALSE
    qc53<-qc53[samples,,drop=FALSE]
}else if(nqc53==2*NS){
    UMI=TRUE
    qc53u<-qc53[paste0("UMI_",samples),,drop=FALSE]
    colnames(qc53u)<-paste0(colnames(qc53),"(UMI)")
    rownames(qc53u)<-sub("^UMI_","",rownames(qc53u))
    qc53<-cbind(qc53[samples,,drop=FALSE],qc53u)
}else{
    stop("star_qc and qc53 do not match")
}

TRIM=dir.exists("fastq_trim")

##collect data for rRNA, globin and phix and fastq_trim/log and duplicates when present
misc<-matrix(NA,NS,5)
colnames(misc)<-c("globin","rRNA","phix","UMI_dup","reads")

for(i in 1:NS){
    SID<-samples[i]
    for(j in 1:3){
        zz<-pipe(paste0("tail -1 ",colnames(misc)[j],"/",SID,".txt |tr -d %"))
        misc[i,j]<-scan(zz,n=1,quiet=TRUE)
        close(zz)
    }
    if(UMI){
        zz<-pipe(paste0("tail -1 star_align/UMI_",SID,"/",SID,"_dup_log.txt |cut -f 6"))
        misc[i,4]<-scan(zz,n=1,quiet=TRUE)*100
        close(zz)
    }
    if(TRIM){
        ##Now it only works for paired
        zz<-pipe(paste0("grep \"Total read pairs processed:\" fastq_trim/log/log.",SID," |cut -f 2 -d: |tr -d ,"))
        misc[i,5]<-scan(zz,n=1,quiet=TRUE)
        close(zz)
    }

}

colnames(star)<-sub(" |","",fixed=TRUE,
                    sub("^\\s+","",perl=TRUE,colnames(star)))

l<-grep("^PCT_",colnames(qc53))[1]
loc<-l:Nqc53
if (UMI){
    loc<-c(loc,loc+Nqc53)
}
qc53<-qc53[,loc,drop=FALSE]

id<-grep("^PCT_",colnames(qc53))
qc53[,id]<-data.matrix(qc53[,id,drop=FALSE])*100

colnames(qc53)[id]<-
    sub(" bases$","",
        paste0("%",gsub("_"," ",
                         gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", perl=TRUE,
                              sub("^PCT_","",colnames(qc53)[id])))))

star<-star[,c(6:7,9:18,25,27,29:31,34),drop=FALSE]
colnames(star)<-c("reads","avg read len","uniq mapped","%uniq mapped","avg mapped len","#splices: Total","#splices: Annotated (sjdb)",
                  "#splices: GT/AG","#splices: GC/AG", "#splices: AT/AC","#splices: Non-canonical", "Mismatch rate per base, %",               
                  "%mapped to multiple loci","%mapped to too many loci" ,"%unmapped: too many mismatches","%unmapped: too short",
                  "%unmapped: other","%chimeric")
star<-data.matrix(sub("%$","",star))
star.info<-dimnames(star)
star<-matrix(as.numeric(star),nrow(star),ncol(star))
dimnames(star)<-star.info

##Now we are working on the columns
colnames(misc)[1:4]<-paste0("%",colnames(misc)[1:4])
qc<-NULL
if(TRIM) qc<-cbind("raw_reads"=misc[,5],"%trim"=star[,1]/misc[,5]*100)
qc<-cbind(qc,misc[,1:3])
if(UMI){
    qc<-cbind(qc,"%UMI_dup"=misc[,4])
}
qc<-cbind(qc,star,qc53)

colnames(qc)[colnames(qc)=="%Mrna"]<-"%mRNA"

write.table(qc,"qc_info.txt",row.names=TRUE,col.names=NA,quote=FALSE,sep="\t")
