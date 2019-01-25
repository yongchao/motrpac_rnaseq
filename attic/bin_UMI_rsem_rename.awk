#!/bin/awk -f
BEGIN{
    FS=OFS="\t"
}
/^@/
!/^@/{
    #Do we need to consider the strand, probably not for rsem, but needed for others
    #remove the not uniquely aligned reads
    umi=gensub("^.*:","","a",$1)
    chr=$3
    pos1=$4
    pos2=$8
    if(and($2,0x001)){
	if (and($2,0xF0c)){
	    next
	}
	if(and(0xfff,compl(and($2,0x002)))){
	    next
	}
	#assuming paired
	if (and($2,0x040)){
	    rname=chr":"pos1"_"pos2":"umi
	}else if(and($2,0x080)){
	    rname=chr":"pos2"_"pos1":"umi	
	}else{
	    print "The mate info is incorrect\n">"/dev/stderr"
	}
    }else{
	#single
	if (and($2,0xF04)){
	    next
	}
	rname=chr":"pos":"umi
    }
    printf rname"\t"
    for(i=2;i<=NF;i++){
	printf "\t"$i
    }
    printf "\n"
}
