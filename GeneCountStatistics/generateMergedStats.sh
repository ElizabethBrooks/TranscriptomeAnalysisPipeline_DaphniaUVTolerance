#!/bin/bash
#Script to generate merged edgeR exact test results and generate top tags
#Usage: bash 
#Usage Ex: bash 

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
outFile=$(basename "$inputsPath"/"$1" | sed 's/\.csv//g')
#Create directory for output files
outDir="$outputsPath"/"$outFile"
outputStats="$outDir"/geneCountStats_cols"$2"to"$3"_"$outFile"
mkdir "$outDir"
mkdir "$outputStats"

#TO DO: Generalize
sed -i 's/"logFC"/"geneID","logFC"/g' topGenesStats_exactTest.csv

cut -f1 -d ',' topGenesStats_exactTest.csv > topGenesTags_exactTest.csv
sed -i 's/"//g' topGenesTags_exactTest.csv
sed -i 's/geneID/Y023/g' topGenesTags_exactTest.csv

paste -d , geneCountStats_cols1to6_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/topGenesTags_exactTest.csv geneCountStats_cols7to12_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/topGenesTags_exactTest.csv geneCountStats_cols13to18_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/topGenesTags_exactTest.csv geneCountStats_cols19to24_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/topGenesTags_exactTest.csv geneCountStats_cols25to30_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/topGenesTags_exactTest.csv geneCountStats_cols31to36_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/topGenesTags_exactTest.csv > mergedTopGenesTags_exactTest.csv
sed -i 's/gene//g' mergedTopGenesTags_exactTest.csv

sed -i 's/"logFC"/"geneID","logFC"/g' stats_exactTest.csv

cut -f1 -d ',' stats_exactTest.csv > statsGenes_exactTest.csv
sed -i 's/"//g' statsGenes_exactTest.csv
sed -i 's/geneID/Sierra/g' statsGenes_exactTest.csv

paste -d , geneCountStats_cols1to6_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/statsGenes_exactTest.csv geneCountStats_cols7to12_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/statsGenes_exactTest.csv geneCountStats_cols13to18_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/statsGenes_exactTest.csv geneCountStats_cols19to24_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/statsGenes_exactTest.csv geneCountStats_cols25to30_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/statsGenes_exactTest.csv geneCountStats_cols31to36_geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned/statsGenes_exactTest.csv > mergedGenes_exactTest.csv
sed -i 's/gene//g' mergedGenes_exactTest.csv