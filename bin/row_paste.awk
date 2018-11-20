#!/bin/awk -f
##row_paste.awk infoid=1 colid=0 <samples.txt
##skip=1, the number of lines to skip
##samples.txt is table delimited file with two columns, one is the the file name and the other is the sample_ID
#if it has only one column, then the file name is also the sample.id
##collect all of the data from many samples
##Assuming one line for the file header
##And Tab delimited file
BEGIN{
    FS=OFS="\t"
    #tmpdir=ENVIRON["MOTRPAC_TMP"]
    tmpdir="/sc/orga/scratch/"ENVIRON["USER"]"/tmp"
    system("mkdir -p "tmpdir)
    outfile=gettmpfilename(tmpdir)
    out1file=gettmpfilename(tmpdir)
    infofile=gettmpfilename(tmpdir)
    info1file=gettmpfilename(tmpdir)
    tmpfile=gettmpfilename(tmpdir)
    headerfile=gettmpfilename(tmpdir)
    
    #print "outfile="outfile"\n" >"/dev/stderr"
    #print "out1file="out1file"\n" >"/dev/stderr"
    #print "infofile="infofile"\n" >"/dev/stderr"
    #print "info1file="info1file"\n" >"/dev/stderr"
    #print "tmpfile="tmpfile"\n" >"/dev/stderr"
    #print "headerfile="headerfile"\n" >"/dev/stderr"
}
NR==1{
    tailh="tail -n +"skip+1" "#with the head line
    tailm="tail -n +"skip+2" "#just the main content

    tailh $1 | getline header
    print header >headerfile
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

    close(headerfile)
    cmd="cut -f "infoid" "headerfile
    #print cmd>"/dev/stderr"
    cmd|getline outhead
    
    sample=$1
    if(NF>=2){
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
