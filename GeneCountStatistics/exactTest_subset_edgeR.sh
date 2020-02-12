#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash exactTest_subset_edgeR.sh countsFolder sample
#Usage Ex: bash exactTest_subset_edgeR.sh GeneCountsAnalyzed_countedCoordinate_htseqHisat2_run1_fullset_run1 R2

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
inFile=$(echo "$inputsPath"/GeneCounts_Merged/"$1"/geneCounts_merged_*_fullset.txt)
noExt=$(echo $inFile | sed 's/\.txt//g')
newFile=$(echo "$noExt".csv)
#Create directory for output files
outDir="$outputsPath"/"$1"
outputStats=$(echo "$outDir"/geneCountStats_"$2")
mkdir "$outDir"
mkdir "$outputStats"
#Convert TXT formatted counts to CSV
sed -e 's/\s\+/,/g' "$inFile" > "$newFile"
#Retrieve selected sample beginning column number
head -1 "$inFile" | tr "\t" "\n" | grep -n "$2" > tmpColNum.txt
colNumStart=$(($(head -1 tmpColNum.txt | cut -d ':' -f1)-1))
colNumEnd=$(($colNumStart+5))
rm tmpColNum.txt
#Perform DE analysis using edgeR and output analysis results to a txt file
Rscript exactTest_edgeR.r "$newFile" $colNumStart $colNumEnd > "$outputStats"/analysisResults.txt
#Rename and move produced filtered table
mv stats_normalizedCounts.csv "$outputStats"/stats_normalizedCounts.csv
mv stats_exactTest.csv "$outputStats"/stats_exactTest.csv
#Rename and move produced plots
mv plotBarsBefore.jpg "$outputStats"/plotBarsBefore.jpg
mv plotMDSBefore.jpg "$outputStats"/plotMDSBefore.jpg
mv plotHeatMapBefore.jpg "$outputStats"/plotHeatMapBefore.jpg
mv plotBarsAfter.jpg "$outputStats"/plotBarsAfter.jpg
mv plotMDSAfter.jpg "$outputStats"/plotMDSAfter.jpg
mv plotHeatMapAfter.jpg "$outputStats"/plotHeatMapAfter.jpg
mv plotBCV.jpg "$outputStats"/plotBCV.jpg
mv plotMD.jpg "$outputStats"/plotMD.jpg
mv plotMA.jpg "$outputStats"/plotMA.jpg
#Move to current outputs folder
cd "$outputStats"
#Make table of GO data for the top tags from exact tests
head -11 stats_exactTest.csv > topGenesStats_exactTest.csv
sed -i 's/"logFC"/"geneID","logFC"/g' topGenesStats_exactTest.csv
cut -f1 -d ',' topGenesStats_exactTest.csv > tmp.csv
sed -i 's/"//g' tmp.csv
head -1 "$ontologyPath" > topGenesGO_exactTest.csv
while IFS= read -r line; do grep $line "$ontologyPath" >> topGenesGO_exactTest.csv; done < "tmp.csv"
rm tmp.csv