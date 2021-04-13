#!/bin/bash

#Script to merge list of gene IDs and pvalues with RBH gene IDs

#RBH file
rbhFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Dmel_Svetec_2016_uvResponsiveGenes_blastp_PA42_v4.1.csv"
#PValue file
pvalueFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Dmel_Svetec_2016_uvResponsiveGenes_pvalues.csv"

#Remove headers
tail -n +2 $queryFileRBH > tmp1.txt
tail -n +2 $dbFileRBH > tmp2.txt

#Pre-clean up
echo "queryHit,dbHit" > $outFileRBH

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q ",$f2" tmp2.txt; then #consensus RBH
		echo "$f1,$f2" >> $outFileRBH
	fi
done < tmp1.txt