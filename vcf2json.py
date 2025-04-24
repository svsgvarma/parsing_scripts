#!/usr/bin/env python3

import os
import sys
import argparse
import json
import re
from pathlib import Path
from collections import defaultdict
from pprint import pprint

def create_vep_dict(info):

    vep_dict = dict()
    vep_list = list()
    if info:
        tmp = info.split('Format:')

        # remove unwanted characters
        rawfields = re.sub('(">)', "", tmp[1])
        fields = rawfields.split('|')
        i = 0
        for field in fields:
            vep_list.append(field)
            vep_dict[field] = i
            i += 1

    return vep_dict, vep_list

def convert_vcf_2_json(vcf):

    out_json = vcf.replace(".vcf", ".json")

    vcf_dict =  defaultdict(dict)
    vcf_dict['header'] =  defaultdict(list)
    vcf_dict['variants'] =  defaultdict(list)

    ##fileformat=VCFv4.2
    ##fileDate=20200418
    ##source=Mutect2,freeBayes
    ##reference=/home/bdelolmo/Desktop/PIPELINE/BUNDLE/ucsc.hg19.fasta
    n_var = 0
    with open (vcf) as f:
        for line in f:
            line = line.rstrip('\n')
            # Parsing header
            if line.startswith('#'):
               # print (line)
                tmp = line.split('=')
                if re.search("##fileformat", line):
                    vcf_dict['header']['fileFormat'] = tmp[1]
                if re.search("##fileDate", line):
                    vcf_dict['header']['fileDate'] = tmp[1]
                if re.search("##source", line):
                    vcf_dict['header']['source'] = tmp[1]
                if re.search("##reference", line):
                    vcf_dict['header']['reference'] = tmp[1]
                if re.search("##contig", line):
                    ###contig=<ID=chr10,length=135534747>
                    line = line.replace("##contig=<ID=", "")
                    line = line.replace("length=", "")
                    line = line.replace(">","")
                    tmp = line.split(",")
                    contig_dict = defaultdict(dict)
                    contig_dict['ID'] = tmp[0]
                    contig_dict['length'] = tmp[1]
                    vcf_dict['header']['contig'].append(contig_dict)

                if re.search("ID=CSQ", line):
                    vep_dict, vep_list = create_vep_dict(line)
                if re.search("ID=CIVIC", line):
                    civic_dict, civic_list = create_vep_dict(line)
                if re.search('##FORMAT', line):
                    line = line.replace("##FORMAT=<ID=", "")
                    line = line.replace(">","")
                    tmp = line.split(',')
                    format_dict = defaultdict(dict)
                    format_id = ''
                    n = 0
                    for item in tmp:
                        tmp2= item.split('=')
                        if n == 0:
                            format_id = item
                        if len(tmp2) >1:
                            format_dict[format_id][tmp2[0]] = tmp2[1]
                        n=n+1
                    vcf_dict['header']['FORMAT'].append(format_dict)
                continue
            else:
                tmp = line.split('\t')
                chr = tmp[0]
                pos = tmp[1]
                id  = tmp[2]
                ref = tmp[3]
                alt = tmp[4]
                qual= tmp[5]
                filter = tmp[6]
                info = tmp[7]
                format_tag = tmp[8]
                format = tmp[9]

                var_name = "var_" + str(n_var)
                vcf_dict['variants'][var_name] = defaultdict(dict)
                vcf_dict['variants'][var_name]['CHROM'] = chr
                vcf_dict['variants'][var_name]['POS'] = pos
                vcf_dict['variants'][var_name]['ID'] = id
                vcf_dict['variants'][var_name]['REF'] = ref
                vcf_dict['variants'][var_name]['ALT'] = alt
                vcf_dict['variants'][var_name]['QUAL'] = qual
                vcf_dict['variants'][var_name]['FILTER'] = filter

                info_list = info.split(';')
                vcf_dict['variants'][var_name]['INFO'] = defaultdict(dict)
                for item in info_list:
                    tmp_item = item.split('=')
                    if item.startswith('CSQ'):
                        vcf_dict['variants'][var_name]['INFO']['CSQ'] = defaultdict(dict)
                        tmp_multidim = item.split(',')
                        n = 0
                        for subitem in tmp_multidim:
                            subitem = subitem.replace("CSQ=", "")
                            conseq_name = "consequence_"+str(n)
                            tmp_subfield = subitem.split('|')
                            num_field = 0
                            for subfield in tmp_subfield:
                                if subfield == "":
                                    subfield = '.'
                                subitem_name = vep_list[num_field]
                                #print (subfield + " " + subitem_name + " " + str(num_field))
                                vcf_dict['variants'][var_name]['INFO']['CSQ'][subitem_name] = subfield
                                num_field=num_field+1
                            n =n+1
                    elif item.startswith('CIVIC'):
                        tmp_multidim = item.split(',')
                        for subitem in tmp_multidim:
                            civic_evidence = "Evidence_"+str(n)
                            tmp_subfield = subitem.split('|')
                            num_field = 0
                            for subfield in tmp_subfield:
                                if subfield == "":
                                    subfield = '.'
                                subitem_name = civic_list[num_field]
                                vcf_dict['variants'][var_name]['INFO']['CIVIC_EVIDENCE'][subitem_name] = subfield
                                num_field=num_field+1
                            n =n+1
                    else:
                        vcf_dict['variants'][var_name]['INFO'][tmp_item[0]] = tmp_item[1]
 
                tmp_format_tag = format_tag.split(':')
                tmp_format = format.split(':')
                j = 0
                for val in tmp_format:
                    vcf_dict['variants'][var_name][tmp_format_tag[j]] = val
                    j=j+1

                n_var = n_var+1
                #chr1	65301110	.	C	T	.	.

    with open(out_json, 'w') as fp:
        json.dump(vcf_dict, fp)


parser = argparse.ArgumentParser(description='Description: Convert VCF to JSON')
parser.add_argument( '--vcf', dest='vcf_file', type = str, 
    help="Input VCF file", required=True)

# Here,obtain arguments
args = parser.parse_args()
vcf = args.vcf_file

convert_vcf_2_json(vcf)

def get_header_from_vcf(vcf):
    header = []
    with open(vcf) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
    return header

