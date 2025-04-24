#!/usr/bin/env python2.7


#python clean_up_input_fasta.py /Volumes/Seagate/Backup-MAC_HD2/proj_Alphonsine/Scripts_blast-positions/Blast-new-genome-pos/Genotype_F2-SNPs-raw_probe.fa /Volumes/Seagate/Backup-MAC_HD2/proj_Alphonsine/Scripts_blast-positions/Blast-new-genome-pos/Genotype_F2-SNPs-raw_probe-clean.fa

import sys
import Bio

from Bio import SeqIO

input_fasta = sys.argv[1]
outputfl = open(sys.argv[2],"w")
for record in SeqIO.parse(input_fasta, "fasta"):
	#hdr=str(" [organism=Lepidium campestre]​")
    fasta_entry = (">{}\n{}\n".format(record.id,record.seq))
    outputfl.write(str(fasta_entry))
outputfl.close()
