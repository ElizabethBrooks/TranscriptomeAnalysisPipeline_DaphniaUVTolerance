#!/bin/bash
#Script to replace gene IDs with entrez IDs with in a raw counts file
#Usage: bash replaceIDs_geneToEntrez.sh rawGeneCounts
#Usage Ex: bash replaceIDs_geneToEntrez.sh cleaned.csv

#Map entrez IDs to gene IDs
tail -n+2 /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezIDs_PA42.csv > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezIDs_tail_PA42.csv
tail -n+2 /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/geneToUniprotIDs_cleaned_PA42.csv > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/geneToUniprotIDs_tail_PA42.csv
sort -t, -k1 /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezIDs_tail_PA42.csv > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezIDs_sorted_PA42.csv
sort -t, -k2 /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/geneToUniprotIDs_tail_PA42.csv > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/geneToUniprotIDs_sorted_PA42.csv
join -1 1 -2 2 -t , -a 1 -e 0 -o '0,1.2,2.1' /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezIDs_sorted_PA42.csv /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/geneToUniprotIDs_sorted_PA42.csv  > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezToGeneIds_PA42.csv
sed '/^gene/d' /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/cleaned_fullSet.csv > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/daphniaRawCounts_onlyEntrez.csv
head -1 /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/cleaned_fullSet.csv > /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/daphniaRawCounts_onlyEntrez_noDups.csv
sort -t ',' -k 1,1 -u /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/daphniaRawCounts_onlyEntrez.csv >> /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/daphniaRawCounts_onlyEntrez_noDups.csv

#Replace gene IDs with entrez IDs
while IFS=, read -r f1 f2 f3
do 
	sed -i "s/$f3,/$f2,/g" /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/cleaned_fullSet.csv
done < /home/mae/Documents/RNASeq_Workshop_ND/PA42_Annotations/uniprotToEntrezToGeneIds_PA42.csv