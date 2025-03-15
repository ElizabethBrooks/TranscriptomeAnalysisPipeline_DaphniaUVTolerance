#!/bin/bash

#Script to merge list of gene IDs and pvalues with RBH gene IDs

#RBH file
#rbhFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Dmel_Svetec_2016_uvResponsiveGenes_blastp_PA42_v4.1.csv"
rbhFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/DEGs_UVBvsCntrl_Tcast_blastp_PA42_v4.1.csv"
#PValue file
#pvalueFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Dmel_Svetec_2016_uvResponsiveGenes_pvalues.csv"
pvalueFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/DEGs_UVBvsCntrl_Tcast_pvalues.csv"

#Output file
#outFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Dmel_Svetec_2016_genePValues.csv"
#outUniqFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Dmel_Svetec_2016_genePValues_uniq.csv"
outFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Tcast_Guo_2019_genePValues.csv"
outUniqFile="/home/mae/Documents/RNASeq_Workshop_ND/CandidateGeneLists/Tcast_Guo_2019_genePValues_uniq.csv"

#Output file
#Remove headers
tail -n +2 $rbhFile > tmp1.csv
tail -n +2 $pvalueFile > tmp2.csv

#Pre-clean up
echo "geneOG,gene,PValue" > $outFile
echo "geneOG,PValue" > $outUniqFile
#sed -i "s/FBpp/FBgn/g" tmp1.csv

#Loop over first set of annotations
while IFS=, read -r f1 f2; do
	#Match gene IDs with PValues
	if grep -q $f1"," tmp1.csv; then
		line=$(grep $f1"," tmp1.csv | head -1)
		echo $line","$f2 >> $outFile
	else
		echo $f1","$f2 >> $outUniqFile
	fi
done < tmp2.csv

#Clean up
rm tmp*.csv