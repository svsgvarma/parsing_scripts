#!/usr/bin/env python3

#python reverse-sequence.py inseq.fa outseq.fa

import sys
import Bio

from Bio import SeqIO

input_fasta = sys.argv[1]
outputfl = open(sys.argv[2],"w")
for record in SeqIO.parse(input_fasta, "fasta"):
	#hdr=str(" [organism=Lepidium campestre]​")
    fasta_entry = (">{}\n{}\n".format(record.id,str(record.seq)[::-1]))
    outputfl.write(str(fasta_entry))
outputfl.close()
