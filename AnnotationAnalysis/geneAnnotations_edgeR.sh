#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash geneAnnotations_edgeR.sh analysisType countsFolder sample
#Usage Ex: bash geneAnnotations_edgeR.sh exactTest genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 R2
#Usage Ex: bash geneAnnotations_edgeR.sh 2WayANOVA_filtered genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 R2

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")

#Create directory for output files
outDir="$inputsPath"/GeneCountsAnalyzed/"$2"/"$3"

#Move to current outputs folder
cd "$outDir"

#Make table of GO data for the top tags from exact tests
cat topTags_"$1".csv > tmpTopTags_"$1".csv
sed -i 's/^logFC/gene,logFC/g' tmpTopTags_"$1".csv
cut -f1 -d ',' tmpTopTags_"$1".csv > tmpTopTagsIds.csv
rm topTagsGO_"$1".csv
#Rertieve top genes GO annotations
while IFS= read -r line; do grep $line "$ontologyPath" >> topTagsGO_"$1".csv; done < "tmpTopTagsIds.csv"
cut -d"," -f2 topTagsGO_"$1".csv | cut -d"^" -f1 > topTagsUniprotIDs_"$1".csv

#Make table of GO data for the top tags from exact tests
cat "$1".csv > tmpGeneIDs_"$1".csv
sed -i 's/^logFC/gene,logFC/g' tmpGeneIDs_"$1".csv
cut -f1 -d ',' tmpGeneIDs_"$1".csv > tmpGeneIDs.csv
cat "$ontologyPath" | cut -d"," -f1,7 | sed "s/\"//g" > tmpGOIDs.csv
rm genesUniprotIDs_"$1".csv 
#Retrieve uniprot IDs from GO annotations
while IFS= read -r line; do grep $line tmpGOIDs.csv >> genesUniprotIDs_"$1".csv; done < "tmpGeneIDs.csv"
cut -d"," -f2 genesUniprotIDs_"$1".csv | cut -d"^" -f1 > uniprotIDs_"$1".csv

#Clean up
rm tmp*.csv