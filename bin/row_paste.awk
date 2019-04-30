#!/bin/awk -f
##row_paste.awk infoid=1 colid=0 <samples.txt
##skip=1, the number of lines to skip
##samples.txt is tab delimited file with two columns, one is the the file name and the other is the sample_ID
#if it has only one column, then the file name is also the sample.id
##collect all of the data from many samples
##Assuming one line for the file header
##skip==-1, means no header
##And Tab delimited file
BEGIN{
    FS=OFS="\t"
    tmpdir="tmpdir" #the currend work folder
    system("mkdir -p "tmpdir)
    outfile=gettmpfilename(tmpdir)
    out1file=gettmpfilename(tmpdir)
    infofile=gettmpfilename(tmpdir)
    info1file=gettmpfilename(tmpdir)
    tmpfile=gettmpfilename(tmpdir)
    headerfile=gettmpfilename(tmpdir)
}
NR==1{
    tailm="tail -n +"skip+2" "#just the main content
    if(skip==-1){
	tailh="head -1 " #a fake head line
    }else{
	tailh="head -"skip+1" |tail -1 "#with the head line
    }

    tailh $1 | getline header
    nf=split(header,headv,"\t")
    if(infoid==0){
	 infoid=1
    }
    if(colid==0){
	colid=nf
    }else if(colid<0){
	colid=nf-colid
    }
    print "infoid="infoid >"/dev/stderr"
    print "colid="colid >"/dev/stderr"

    if(skip>=0){
	outhead=headv[colid]
    }else{
	outhead="field"
    }
    
    sample=$1
    if(NF>=2){#the last column may become the label
	sample=$2
    }

    outhead=outhead"\t"sample
    cmd=tailm $1"|cut -f "infoid" >"infofile
    #print cmd >"/dev/stderr"
    system(cmd) 
    cmd=tailm $1"|cut -f "colid"|paste "infofile" - >"outfile
    #print cmd >"/dev/stderr"
    system(cmd)
}	
NR!=1{
    sample=$1
    if(NF>=2){
	sample=$2
    }
    outhead=outhead"\t"sample
    system(tailm $1"|cut -f "infoid" >"info1file) 
    status=system("cmp "info1file" "info1file">/dev/null")
    if(status!=0){
       print "The file "$1 " has differenct info\n"
       exit 1
   }
   system(tailm $1"|cut -f "colid" >"out1file) 
##sequentially paste rather than parrall paste to make programming easier but less efficient
   system("paste "outfile" "out1file" >"tmpfile)
   system("sleep 1") ##sleep one second
   system("mv "tmpfile" "outfile)
}
END{
    print outhead
    system("cat "outfile)
    #remove the temporary files
    system("rm -f "outfile)
    #print outfile >"/dev/stderr" 
    system("rm -f "out1file)
    system("rm -f "infofile)
    system("rm -f "info1file)
    system("rm -f "tmpfile)
    system("rm -f "headerfile)
}
function gettmpfilename(tmpdir,  res,cmd){
    cmd="mktemp --tmpdir="tmpdir
    cmd|getline res
    close(cmd)
    return res
}
