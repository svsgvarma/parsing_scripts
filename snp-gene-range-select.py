#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Parse polymarker input and output a fasta file for blast.

./snp-gene-range-select.py Asativa_sang.v1.1_pars.bed /mnt/c/proj/REF/gene_annotation/

#cat Asativa_sang.v1.1.gff3 | grep "gene" | awk '{split($9,a,";"); print $1"\t"$4"\t"$5"\t"a[1]}' FS="\t" OFS="\t" | sed 's/ID=//g' > Asativa_sang.v1.1_pars.bed

"""

import numpy as np
import sys
import subprocess

def main():
    
    ##############################################
    #####--Function to search SNP data ---
    def srchgene_ranges_find(GRdict,SNPpos):
        SNPpos_sp = SNPpos.split(":")
        idx_rng=[]
        for keys, values in GRdict.items():
            keys_sp1 = keys.split(",")[0].split("-")
            idx_rng.append(int(SNPpos_sp[1]) in range(int(keys_sp1[0]), int(keys_sp1[1])+1))
        idx_pos = np.where(idx_rng)[0]
        snpinlst = [list(GRdict.items())[i] for i in idx_pos]
        snpinlst_indx = [i for i in snpinlst if(i[1] == SNPpos_sp[0])]
        return snpinlst_indx
    ###############################################
    #### test SNPs
    #GRdict = { "167771950-167772015,gene01" : "chr1", "167772352-167772415,gene1" : "chr1", "167772160-167772223,gene22" : "chr2", "167772288-167772351,gene23" : "chr1", "167772224-167772255,gene2" : "chr3", "167772335-167772355,gene432" : "chr4"}
    #SNPpos = 167772341
    #SNPpos = "chr1:167772341"
    #SNPpos = "chr1:167772355"
    #SNPpos = "chr2:167772054"
    #SNPpos = "chr1:167772054"
    #### Actual SNP 
    SNPpos = "chr7D:223545516"
    #SNPpos = "chr7D:223546409"

    #polymarker_input = sys.argv[1]
    #workdir=sys.argv[2]
   
    workdir="/home/varma/proj/REF/Avena-sativa_Sangv1.1/gene_annotation/"
    polymarker_input = "Asativa_sang.v1.1.gff3"

    GRdict={}
    cmdFls1 = subprocess.check_output("cat "+workdir+polymarker_input+" | grep 'gene' ",shell=True)
    #print(cmdFls1)
    for line in cmdFls1.strip().decode().split("\n"):
        lines = line.split("\t")
        keys_all=str(lines[3]+"-"+lines[4]+","+lines[8].split(";")[0].split("=")[1])
        GRdict[keys_all] = str(lines[0])

    snpinlst_fltr = srchgene_ranges_find(GRdict,SNPpos)
    ############ -------- #############
    snpinlst_fltr_cln=[]; snpinlst_fltr_updowncln=[]
    if (snpinlst_fltr == []):
        for i in range(0, 1000000, 100):
            SNPpos_num = SNPpos.split(":")[1]
            ##--UP
            SNPpos_up = str(SNPpos.split(":")[0])+":"+str(int(SNPpos_num)-i)
            snpinlst_fltr_up = srchgene_ranges_find(GRdict,SNPpos_up)
            ##--DOWN
            SNPpos_down = str(SNPpos.split(":")[0])+":"+str(int(SNPpos_num)+i)
            snpinlst_fltr_down = srchgene_ranges_find(GRdict,SNPpos_down)
            if (snpinlst_fltr_up != []):
                True
                snpinlst_fltr_updowncln.append(str("UP:-")+str(i)+";"+str(snpinlst_fltr_up[0][1]+":"+snpinlst_fltr_up[0][0]))
                break
            elif (snpinlst_fltr_down != []):
                True
                snpinlst_fltr_updowncln.append(str("DOWN:+")+str(i)+";"+str(snpinlst_fltr_down[0][1]+":"+snpinlst_fltr_down[0][0]))
                break
            elif (snpinlst_fltr_down == []):
                True
                snpinlst_fltr_updowncln.append(str("."))
    else:
        snpinlst_fltr_cln.append(str(snpinlst_fltr[0][1]+":"+snpinlst_fltr[0][0]))
    
    if(snpinlst_fltr_cln!=[]):
        print(snpinlst_fltr_cln[0])
    elif(snpinlst_fltr_updowncln!=[] ):
        print(snpinlst_fltr_updowncln[0])

    #######
    print("script done!")
    return 0

if __name__ == '__main__':
    main()

