#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash statistics_exactTest_subset.sh countsFolder startColPos endColPos
#Usage Ex: bash statistics_exactTest_subset.sh GeneCountsAnalyzed_countedCoordinate_htseqHisat2_run1_fullset_run1 1 6

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
inFile="$inputsPath"/GeneCounts_Merged/"$1"/geneCounts_merged_*_fullset.txt
#Create directory for output files
outDir="$outputsPath"/"$1"
outputStats="$outDir"/geneCountStats_cols"$2"to"$3"
mkdir "$outDir"
mkdir "$outputStats"

#Perform DE analysis using edgeR and output analysis results to a txt file
Rscript exactTest_edgeR.r "$inFile" $2 $3 > "$outputStats"/analysisResults.txt
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
head -1 "../gene.Blast2GO.merged.csv" > topGenesGO_exactTest.csv
while IFS= read -r line; do grep $line "../gene.Blast2GO.merged.csv" >> topGenesGO_exactTest.csv; done < "tmp.csv"
rm tmp.csv