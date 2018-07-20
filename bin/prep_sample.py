#!/usr/bin/env python

##How to make sure the python is this one
#/sc/orga/projects/haghif01a/local/mini_conda/python3/bin/python

##input: the fastq files under the fastq (or fastq_raw) folder
#output: samples_info.txt data.frame with columns SID, barcode, lane, pairs, fastqfile, totally data dependent
#barcodes was looking at the random 100 reads and uses the top 2 most frequent one.

import glob
from collections import defaultdict

def freq(x):
    d = defaultdict(int)
    for i in x:
        d[i] += 1
    return(d)

def sort_freq(d):
    res=defaultdict(int)    
    for w in sorted(d, key=d.get, reverse=True):
        res[w]=d[w]
    return res    

def top2_freq(d,info,n,res):
    sd=sort_freq(d)
    #the first entry
    key=list(sd)[0]
    frac=sd[key]/(n+0.0)
    res[info+"1"]="%s:%0.2f" % (key,frac)
    #the second entry
    if len(sd)>1:
        key=list(sd)[1]
        frac=sd[key]/(n+0.0)
    else:
        key=""
        frac=0.0
    res[info+"2"]="%s:%0.2f" % (key,frac)

def fastq_info(x):
    "Obtain the lane and squenece length, read and barcode info"
    "Currently look for the top 100 reads, future plan is to do random 100 reads"
    import gzip
    n=0
    seqid=defaultdict(int)
    lane=defaultdict(int)
    r=defaultdict(int)
    #The above three have to be unique
    bc=defaultdict(int)
    seqlen=defaultdict(int)
    with gzip.open(x,'rt') as f:
        for line in f:
            line=line.strip() #important to remove the white spaces and \n
            n+=1
            if n%4==1:
                #work on the first line to get the eensential info
                L=line.split()
                L1=L[0].split(":")
                lane[L1[-4]]+=1
                seqid[L1[-5]]+=1
                L2=L[1].split(":")
                r[L2[0]]+=1
                bc[L2[3]]+=1
            elif n%4==2:
                seqlen[str(len(line))]+=1
            if n>=400:
                #Find the top 100 reads
                break
    #get the most frequent one
    if len(seqid) !=1 or (len(lane) !=1 or len(r) !=1):
        print("Error: one of the reads or sequecing projectid or sequencing lane number is not uique in file "+x)
        exit(1)
    info={}
    info["r"]=list(r)[0]
    info["lane"]=list(lane)[0]
    info["seqid"]=list(seqid)[0]
    #print the top 2 entries
    n=n//4
    top2_freq(bc,"bc",n,info)
    top2_freq(seqlen,"seqlen",n,info)
    return(info)

def prep_samples(xv):
   SR={}
   for x in xv:
       r=x[-2:] #R1 or R2
       s=x[:-3] #Sample name
       if s in SR:
          SR[s]=[r]+SR[s]
       else:
          SR[s]=[r]
   samples=defaultdict(dict)  
   for s,r in SR.items():
       OK=True
       R1=r.count('R1')
       R2=r.count('R2')
       if R1==1:
          if R2==1:
             samples[s]['R']='R1,R2'
          elif R2==0:
             samples[s]['R']='R1'
          else: OK=False
       else:
          OK=False
       if not OK:   
          print("Error: no consistent R1 and R2 files for sample "+s)
   ##obtain the fastq_info
   for s in samples:
       si=samples[s]
       if si['R']=='R1':
           info=fastq_info(s+"_R1.fastq.gz")
           if info["r"]!='1':
               print("Error: "+s+" the read info entry in R1 file is not 1")
               exit(1)
           for i in info:
               si[i]=info[i]
       else:
           info1=fastq_info(s+"_R1.fastq.gz")
           info2=fastq_info(s+"_R2.fastq.gz")
           if info1["r"]!='1' and info2["r"]!='2':
               print("Error: "+s+" the read info entry in R1 or R2 files is not 1 (for R1) or 2 (for R2)")
               exit(1)
           if not info1["lane"]==info2["lane"]:
               print("Error: lane number entry of R1 and R2 files for sample "+s+" is not the same")
               exit(1)
           si["lane"]=info1["lane"]
           
           if not info1["seqid"]==info2["seqid"]:
               print("Error: the seq project id for R1 and R1 files of sample "+s+" is not the same")
               exit(1)
           si["seqid"]=info1["seqid"]
           keys=["bc1","bc2","seqlen1","seqlen2"]
           for j in keys:
               si[j]=info1[j]+","+info2[j]
   return samples

#main
files=glob.glob("*.fastq.gz")
fastqc_samples=[x[:-len('.fastq.gz')] for x in files]
samples=prep_samples(fastqc_samples)
keys=["R","lane","seqid","bc1","bc2","seqlen1","seqlen2"]
print("SID\t"+'\t'.join(keys))
for s in samples:
    si=samples[s]
    sout=s;
    for j in keys:
        sout+="\t"+si[j]
    print(sout)
