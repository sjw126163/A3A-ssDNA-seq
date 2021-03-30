#!/usr/bin/env python3
import sys
import re
import pysam
infile1 =pysam.AlignmentFile(sys.argv[1],"rb")
name = sys.argv[1][:-4]+"_C_T_position.bed"
outfile1 = open(name, "w")
p1 = re.compile(r"\D+")
p2 = re.compile(r"\d+")
for line in infile1:
    header = line.query_name
    mapq = line.mapping_quality
    chr_name = str(line.reference_name)
    start_pos = line.reference_start
    if mapq <= 10:
        continue

    if line.is_reverse:
        strand = "-"
    else:
        strand = "+"
    seq = line.query_alignment_sequence
    if line.has_tag('MD'):
        mismatch = line.get_tag('MD')+"~"
        refBases = p1.findall(mismatch)
        refPositions = p2.findall(mismatch)
        i = 0
        refPos = start_pos
        readPos = 0
        while i < len(refBases)-1:
            refBase = refBases[i]
            pos = int(refPositions[i])
            readPos = readPos + pos
            refPos = refPos + pos

            #print(seq[readPos],refBase)
            if refBase == "C" and seq[readPos] == "T" :
                print(chr_name, refPos, refPos + 1, header, mapq, strand, sep="	", file=outfile1)

            if refBase[0] != "^":
                readPos = readPos + 1
                refPos = refPos + 1
            else:
                refPos = refPos + len(refBase)-1

            i = i+1
infile1.close()
outfile1.close()

        
