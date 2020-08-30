#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash exactTest_geneAnnotations_edgeR.sh countsFolder sample
#Usage Ex: bash exactTest_geneAnnotations_edgeR.sh genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 R2

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
outDir="$inputsPath"/GeneCountsAnalyzed/"$1"/"$2"

#Move to current outputs folder
cd "$outDir"

#Make table of GO data for the top tags from exact tests
head -11 stats_exactTest.csv > topGenesStats_exactTest.csv
sed -i 's/^logFC/gene,logFC/g' topGenesStats_exactTest.csv
cut -f1 -d ',' topGenesStats_exactTest.csv > tmpTopGeneIds.csv
rm topGenesGO_exactTest.csv
#Rertieve top genes GO annotations
while IFS= read -r line; do grep $line "$ontologyPath" >> tmpTopGenesGO_exactTest.csv; done < "tmpTopGeneIds.csv"

#Make table of GO data for the top tags from exact tests
cat stats_exactTest.csv > tmpGeneIDs_exactTest.csv
sed -i 's/^logFC/gene,logFC/g' tmpGeneIDs_exactTest.csv
cut -f1 -d ',' tmpGeneIDs_exactTest.csv > tmpGeneIDs.csv
cat "$ontologyPath" | cut -d"," -f1,7 | sed "s/\"//g" > tmpGOIDs.csv
rm genesUniprotIDs_exactTest.csv 
#Retrieve uniprot IDs from GO annotations
while IFS= read -r line; do grep $line tmpGOIDs.csv >> genesUniprotIDs_exactTest.csv; done < "tmpGeneIDs.csv"
cut -d"," -f2 genesUniprotIDs_exactTest.csv | cut -d"^" -f1 > uniprotIDs_exactTest.csv

#Clean up
rm tmp*.csv