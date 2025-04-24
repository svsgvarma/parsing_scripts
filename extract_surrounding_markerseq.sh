#!/usr/bin/sh

#alpha version 0.1

############################
############################
# Use pars ploymarker flaking sequence script

#########################################################################
# for short range 
# Pangenome path:

INDIR=/mnt/c/proj/proj_extract-marker-seq/data/
OUTDIR=/mnt/c/proj/proj_extract-marker-seq/work/

REF=${INDIR}"genomic.fasta"
INFILE=${INDIR}"selected_makrers.csv"
OUTFILE=${OUTDIR}"selected_makrers_seq-updown50.tsv"
updownrange=50

rm ${OUTFILE}

# Define the header
header="Sno\tChr\tPos\tPos_updown50\tSeq_updown50"

echo -e "$header" > "${OUTFILE}"

while read p; do
  #echo "$p"
  posNO=$(echo $p | awk -F"," '{print $1}' | sed 's/"//g')
  chrID=$(echo $p | awk -F"," '{print $2}' | sed 's/-/./g'| sed 's/"//g')
  markerpos=$(echo $p | awk -F"," '{print $3}' | grep -o '[0-9]*')

  ##############
  #echo ${posNO}
  #echo ${chrID}
  #echo ${markerpos}
  ##############
  # Define position ranges 
  markerup=$(expr ${markerpos} - ${updownrange} )
  markerdown=$(expr ${markerpos} + ${updownrange})
  markerposrange="${chrID}:${markerup}-${markerdown}"
  #echo ${markerposrange}

  ##############
  # Search in the reference genome 
  # Use sed to replace newlines with empty string except lines starting with '>'
  findstring=$(samtools faidx ${REF} ${markerposrange} | sed -r 's/>/>'${markerpos}'_@_/g')

  # Use awk to Process the FASTA format with tab separation after headers and join all the sequences
  findstring_output=$(echo "$findstring" | awk '/^>/ {if (NR > 1) printf("\t"); printf("%s\t", $0); next} {printf("%s", $0)} END {printf("\n")}' | tr '[:lower:]' '[:upper:]')

  # Print the result
  echo -e "$posNO\t$chrID\t$markerpos\t$findstring_output" >> ${OUTFILE}

done < ${INFILE}

#########################################################################
