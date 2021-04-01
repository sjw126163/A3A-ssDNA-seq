#!/usr/bin/env python3
import sys
import re
import pysam
infile1 =pysam.AlignmentFile(sys.argv[1],"rb")
name = sys.argv[1][:-4]+"_C_T_position_sjw.bed"
outfile1 = open(name, "w")
p1 = re.compile(r"\D+")
p2 = re.compile(r"\d+")
for line in infile1:
    header = line.query_name
    mapq = line.mapping_quality
    chr_name = str(line.reference_name)
    start_pos = line.reference_start
    seq = line.query_alignment_sequence
    if mapq <= 10:
        continue
    if line.is_read1:
    	read_info ="R1"
    	
    	if line.has_tag('MD'):
    		mismatch = line.get_tag('MD')+"~"
    		refBases = p1.findall(mismatch)
    		refPositions = p2.findall(mismatch)
    		i =0
    		refPos = start_pos
    		readPos = 0
    		while i < len(refBases)-1:
    			refBase = refBases[i]
    			pos = int(refPositions[i])
    			readPos = readPos + pos
    			refPos = refPos + pos
    			if refBase == "G" and seq[readPos] == "A" and line.is_reverse:
    				print(chr_name, refPos, refPos + 1, header, mapq, "-", read_info,sep="\t", file=outfile1)
    			if refBase == "C" and seq[readPos] == "T" and not line.is_reverse:
    				print(chr_name, refPos, refPos + 1, header, mapq, "+", read_info,sep="\t", file=outfile1)
    			if refBase[0] != "^":
    				readPos = readPos + 1
    				refPos = refPos + 1
    			else:
    				refPos = refPos + len(refBase)-1
    			i =i+1
                
    if line.is_read2:
    	read_info ="R2"
    	if line.has_tag('MD'):
    		mismatch = line.get_tag('MD')+"~"
    		refBases = p1.findall(mismatch)
    		refPositions = p2.findall(mismatch)
    		i =0
    		refPos = start_pos
    		readPos = 0
    		while i < len(refBases)-1:
    			refBase = refBases[i]
    			pos = int(refPositions[i])
    			readPos = readPos + pos
    			refPos = refPos + pos
    			if refBase == "C" and seq[readPos] == "T" and line.is_reverse:
    				print(chr_name, refPos, refPos + 1, header, mapq, "+", read_info,sep="\t", file=outfile1)
    			if refBase == "G" and seq[readPos] == "A" and not line.is_reverse:
    				print(chr_name, refPos, refPos + 1, header, mapq, "-", read_info,sep="\t", file=outfile1)
    			if refBase[0] != "^":
    				readPos = readPos + 1
    				refPos = refPos + 1
    			else:
    				refPos = refPos + len(refBase)-1
    			i =i+1
    	
    	
infile1.close()
outfile1.close()

        
