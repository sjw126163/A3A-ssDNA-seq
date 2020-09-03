#!/usr/bin/env python3
import sys
import re
import pysam
infile1 =pysam.AlignmentFile(sys.argv[1],"rb")
name =sys.argv[1][:-4]+"_C_T_position.bed"
outfile1 =open(name,"w")
for line in infile1:
    if line.is_reverse ==True:
    	mappos = "-"
    if line.is_reverse ==False:
    	mappos ="+"
    header =line.query_name
    mapq =line.mapping_quality
    chr_name =str(line.reference_name)
    start_pos =line.reference_start
    seq =line.seq
    length =0
    storage =""
    if line.has_tag('MD'):
        mismatch =line.get_tag('MD')+"~"
        dic1 ={"0","1","2","3","4","5","6","7","8","9"}
        dic2 ={"A","N","T","G","a","t","g"}
        for character in mismatch:
            if character in dic1:
                storage += character
            if character not in dic1 and storage !="":
                length += int(storage)
            if character =="^":
                storage =""
                regexDe =re.search(r'\^+[A-Z]{1,}',mismatch)
                if regexDe:
                    count =int(regexDe.end(0)-regexDe.start(0)-1)
                    seq =seq[0:length]+str("0"*count)+seq[length:]
                    end =regexDe.end(0)
                    mismatch =mismatch[end:]
            if character in dic2:
                length +=1
                storage =""
            if character =="C":
                storage =""
                if seq[length] =="T":
                    pos1 =start_pos+length
                    pos2 =pos1+1
                    print(chr_name,pos1,pos2,header,mapq,mappos,sep="	",file =outfile1)
                length +=1
            if character =="c":
                storage =""
                if seq[length] =="t":
                    pos1 =start_pos+length
                    pos2 =pos1+1
                    print(chr_name,pos1,pos2,header,mapq,mappos,sep="	",file =outfile1)
                length +=1
infile1.close()
outfile1.close()

        