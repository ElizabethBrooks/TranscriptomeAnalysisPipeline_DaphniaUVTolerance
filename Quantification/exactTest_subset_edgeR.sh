#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash exactTest_subset_edgeR.sh countsFolder sample
#Usage Ex: bash exactTest_subset_edgeR.sh genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 R2

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
inFile=$(echo "$inputsPath"/GeneCountsAnalyzed/"$1"/cleaned.csv)
noExt=$(echo $inFile | sed 's/\.txt//g')
newFile=$(echo "$noExt".csv)

#Create directory for output files
outDir="$inputsPath"/GeneCountsAnalyzed/"$1"/"$2"
mkdir $outDir

#Convert TXT formatted counts to CSV
#sed -e 's/\s\+/,/g' "$inFile" > "$newFile"
cat "$inFile" > "$newFile"
#Retrieve selected sample beginning column number
#head -1 "$inFile" | tr "\t" "\n" | grep -n "$2" > tmpColNum.txt
head -1 "$newFile" | tr "," "\n" | grep -n "$2" > tmpColNum.txt
colNumStart=$(($(head -1 tmpColNum.txt | cut -d ':' -f1)-1))
colNumEnd=$(($colNumStart+5))

#Clean up
rm tmpColNum.txt

#Perform DE analysis using edgeR and output analysis results to a txt file
Rscript exactTest_edgeR.r "$newFile" $colNumStart $colNumEnd > "$outDir"/analysisResults.txt

#Rename and move produced normalized counts table and exact test stats
mv stats_normalizedCounts.csv "$outDir"/stats_normalizedCounts.csv
mv stats_exactTest.csv "$outDir"/stats_exactTest.csv
#Rename and move produced plots
mv plotBarsBefore.jpg "$outDir"/plotBarsBefore.jpg
mv plotMDSBefore.jpg "$outDir"/plotMDSBefore.jpg
mv plotHeatMapBefore.jpg "$outDir"/plotHeatMapBefore.jpg
mv plotBarsAfter.jpg "$outDir"/plotBarsAfter.jpg
mv plotMDSAfter.jpg "$outDir"/plotMDSAfter.jpg
mv plotHeatMapAfter.jpg "$outDir"/plotHeatMapAfter.jpg
mv plotBCV.jpg "$outDir"/plotBCV.jpg
mv plotMD.jpg "$outDir"/plotMD.jpg
mv plotMA.jpg "$outDir"/plotMA.jpg

#Move to current outputs folder
#cd "$outDir"

#Fix formatting and headers for the normalized counts table and exact test stats
#sed -i 's/"//g' stats_normalizedCounts.csv
#sed -i 's/"//g' stats_exactTest.csv
#sed -i -e 's/^$2_VIS_Pool1/gene,$2_VIS_Pool1/' stats_normalizedCounts.csv
#sed -i -e 's/^logFC/gene,logFC/' stats_exactTest.csv

#Make table of GO data for the top tags from exact tests
#head -11 stats_exactTest.csv > topGenesStats_exactTest.csv
#sed -i 's/^logFC/gene,logFC/g' topGenesStats_exactTest.csv
#cut -f1 -d ',' topGenesStats_exactTest.csv > tmpTopGeneIds.csv
#head -1 "$ontologyPath" > topGenesGO_exactTest.csv
#Rertieve top genes GO annotations
#while IFS= read -r line; do grep $line "$ontologyPath" >> topGenesGO_exactTest.csv; done < "tmpTopGeneIds.csv"

#Make table of GO data for the top tags from exact tests
#cut -f1 -d ',' stats_exactTest.csv > tmpGeneIDs_exactTest.csv
#cat "$ontologyPath" | cut -d"," -f1,7 | sed "s/\"//g" > tmpGOIDs.csv
#head -1 "$ontologyPath" | cut -d"," -f1,7 | sed "s/\"//g" > genesUniprotIDs_exactTest.csv 
#Retrieve uniprot IDs from GO annotations
#while IFS= read -r line; do grep $line tmpGOIDs.csv >> genesUniprotIDs_exactTest.csv; done < "tmpGeneIDs_exactTest.csv"

#Clean up
#rm tmp*.csv