#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla
##work ont the qc53.plot
#the folder should be in the parent folder of qc53, the project folder
options(stringsAsFactors=FALSE)
.DEBUG=interactive()
if(.DEBUG){
    Args<-"sample_info.txt"
}else{
    Args <-commandArgs(TRUE)
}

if(length(Args)>=1){
    file.samples<-Args[1]
}else{
    file.samples<-"sample_info.txt"
}

if(length(Args)>=2){
    proj<-Args[2]
}else{
    proj<-"qc53"
}


sample.info<-read.delim(file.samples,sep="\t",head=TRUE,row.names=1,check.names=FALSE)
samples<-row.names(sample.info)
if ("label" %in% colnames(sample.info)){
    labels<-sample.info[,"label"]
    names(labels)<-samples
}else{
    labels<-samples
}
    
n<-length(samples)
x<-NULL
head<-NULL
for(i in 1:length(samples)){
    filei<-file.path("qc53",paste0(samples[i],"_0-0k.RNA_Metrics"))
    x<-cbind(x,scan(filei,"",sep="\t",skip=7,nline=1,quiet=TRUE))
    headi<-scan(filei,"",sep="\t",skip=6,nline=1,quiet=TRUE)
    if(is.null(head)){
        head<-headi
    }else{
        if(sum(head!=headi)>0){
            stop("Error:",samples[i])
        }
    }
}
dimnames(x)<-list(head,samples)
write.table(x, file.path("qc53",paste0(proj,".txt")),
            sep="\t",col.names=NA,row.names=TRUE)
            
types<-paste0(c("0-0","0-1","1-2","2-3","3-4","4-6","6-10","10-0"),"k")
if(n<=8){
    cols<-palette()[1:n] #only eight colors
}else{
    #if(n<=14){
    #library(RColorBrewer)
    #cols<-c("black","gray",brewer.pal(n-2,"Paired"))
#}else{
    pal<-matrix(c(0,0,0,
                  1,0,0,
                  0,1,0,
                  0,0,1,
                  1,1,0,
                  0,1,1,
                  1,0,1),ncol=3,byrow=TRUE)
    seg<-1:ceiling(n/nrow(pal))
    cols<-matrix("",nrow=length(seg),ncol=7)
    for(i in seg){
        cols[i,]<-rgb(pal[,1],pal[,2],pal[,3],i/length(seg))
    }
    cols<-c(cols)[1:n]
}
pdf(file.path("qc53",paste0(proj,".pdf")),
    width=8.5,height=11)
par(mfrow=c(2,1))
for(j in 1:length(types)){
    cat("Running", types[j],"\n")
    x<-NULL
    samplesj<-NULL
    for(i in 1:length(samples)){
        filei<-file.path("qc53",paste0(samples[i],"_",types[j],".RNA_Metrics"))
        zz<-pipe(paste0("wc -l ",filei))
        nl<-scan(zz,n=1,quiet=TRUE)
        close(zz)
        if(nl==10){
            ##no data, we need to skip
            next
        }
        xi<-read.delim(filei,sep="\t",skip=10,head=TRUE,nrows=101,row.names=1)
        if(length(samplesj)==0){
            x<-xi
        }else{
            if(sum(rownames(xi)!=rownames(x))!=0){
                stop("Error in reading file",filei,"\n")
            }
            x<-cbind(x,xi)
        }
        samplesj<-c(samplesj,samples[i])
    }

    if(is.null(x)){
        next #with no data for any sample, skip
    }
    colnames(x)<-samplesj
    x<-data.matrix(x)
    pos<-as.numeric(rownames(x))
    xrng<-range(pos)
    yrng<-range(x)
    main<-types[j]
    if(main=="0-0k"){
        main<-"default (all)"
    }else if(main=="10-0k"){
        main<-"10k+"
    }
    
    plot(xrng,yrng,type="n",xlab="5prime to 3prime positions",ylab="normalized average read coverage",main=main)
    for(i in 1:ncol(x)){
        lines(pos,x[,i],col=cols[i])
    }
    legend("bottom",
           legend=labels[samplesj],col=cols,lty=1,cex=0.6)
}

dev.off()
