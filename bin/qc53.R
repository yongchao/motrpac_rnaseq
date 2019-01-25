#!/usr/bin/env Rscript
##Note --vanilla is supported only for the actual path
##Input: qc53/{$SID}.RNA_Metrics" and possibly qc53/UMI_{$SID}.RNA_Metrics
##output: qc53.plot and qc53.txt
                                        
##the program should be run in the project folder, the parent folder of qc53

options(stringsAsFactors=FALSE)
.DEBUG=interactive()
if(.DEBUG){
    Args<-"samples_all"
}else{
    Args <-commandArgs(TRUE)
}

if(length(Args)>=1){
    file.samples<-Args[1]
}else{
    file.samples<-"samples" #samples_all
}

if(length(Args)>=2){
    proj<-Args[2]
}else{
    proj<-"qc53"
}


samples<-scan(file.samples,"",quiet=TRUE)
files<-file.path("qc53",paste0(samples,".RNA_Metrics"))
n<-length(samples)
x<-NULL
head<-NULL
for(i in 1:length(samples)){
    x<-cbind(x,scan(files[i],"",sep="\t",skip=7,nline=1,quiet=TRUE))
    headi<-scan(files[i],"",sep="\t",skip=6,nline=1,quiet=TRUE)
    if(is.null(head)){
        head<-headi
    }else{
        if(sum(head!=headi)>0){
            stop("Error:",samples[i])
        }
    }
}
dimnames(x)<-list(head,samples)
write.table(x, paste0(proj,".txt"),
            sep="\t",col.names=NA,row.names=TRUE)
            
##Plot every 8 samples, with the 
cols<-palette()[1:8] #only eight colors
pdf(paste0(proj,".pdf"),
    width=8.5,height=11)
par(mfrow=c(2,1))
nP<-ceiling(length(samples)/8)
for(i in 1:nP){
    locK<-((i-1)*8+1):min(i*8,length(samples))
    ##checking no data case
    x<-NULL
    for(k in locK){
        zz<-pipe(paste0("wc -l ",files[k]))
        nl<-scan(zz,n=1,quiet=TRUE)
        close(zz)
        if(nl==10){
        ##no data, we need to skip
            next
        }
        xk<-read.delim(files[k],sep="\t",skip=10,head=TRUE,nrows=101,row.names=1)
        if(is.null(x)){
            x<-xk
        }else if(sum(rownames(x)!=rownames(xk))!=0){
            stop("Error in reading file",files[k],"\n")
        }else{    
            x<-cbind(x,xk)
        }
        colnames(x)[ncol(x)]<-samples[k]
    }
    if(!is.null(x)){
        x<-data.matrix(x)
        pos<-as.numeric(rownames(x))
        xrng<-range(pos)
        yrng<-range(x)
        plot(xrng,yrng,type="n",xlab="5prime to 3prime positions",ylab="normalized average read coverage")
        for(i in 1:ncol(x)){
            lines(pos,x[,i],col=cols[i])
        }
        legend("bottom",
               legend=colnames(x),col=cols,lty=1,cex=0.6)
    }
}
dev.off()
