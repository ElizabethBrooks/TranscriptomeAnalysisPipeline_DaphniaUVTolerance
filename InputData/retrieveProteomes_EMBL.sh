#!/bin/bash
#Script to retrieve proteome fasta files from EMBL-EBI
# based on input tab delimited proteome ID file downloaded
# from an EMBL-EBI reference proteome search
#Usage: bash retrieveProteomes_EMBL.sh tabDelimitedIDFile
#Usage ex: bash retrieveProteomes_EMBL.sh proteomes-drosophila-filtered-reference_yes.tab

#Retrieve input proteome ID list
cut -f 1 "$1" > proteomeIDs.txt
sed -i 1d proteomeIDs.txt #Remove header

while read pID; do
  echo "$pID"
  curl https://www.uniprot.org/uniprot/?query=proteome:"$pID"&format=fasta
  #curl -X GET --header 'Accept:application/json' 'https://www.ebi.ac.uk/proteins/api/proteomes/"$pID"'
done < proteomeIDs.txt
