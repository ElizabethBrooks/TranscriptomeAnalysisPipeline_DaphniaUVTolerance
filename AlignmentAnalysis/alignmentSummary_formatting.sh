#!/bin/bash
#Script to prepare alignment summary matrices
#Usage: bash alignmentSummary_formatting.sh analysisMethod runNum
#Usage Ex: bash alignmentSummary_formatting.sh tophat2 3

#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory and file
outDir="$outputsPath"/AlignmentsAnalyzed
outFile="$outDir"/alignmentSummarized_"$1"_run"$2"_formatted.csv
outFileTmp="$outDir"/alignmentSummarized_"$1"_run"$2"_tmp.csv
#Set input path
inFile="$outDir"/alignmentSummarized_"$1"_run"$2".csv
#Remove percent symbols
sed "s/%//g" "$inFile" > "$outFile"
#Remove excess tagging
sed -i "s/140327_I481_FCC3P1PACXX_L2_//g" "$outFile"
sed -i "s/140327_I481_FCC3P1PACXX_L3_//g" "$outFile"
sed -i "s/140327_I481_FCC3P1PACXX_L4_//g" "$outFile"
#Re-format pool tag
sed -i "s/Pool_1/Pool1/g" "$outFile"
sed -i "s/Pool_2/Pool2/g" "$outFile"
sed -i "s/Pool_3/Pool3/g" "$outFile"
#Add matrix header
head -1 "$outFile" > "$outFileTmp"
#Re-order UV samples
grep "Pool1.*UV" "$outFile" | sed "s/UV/UV_Pool1/g" >> "$outFileTmp"
grep "Pool2.*UV" "$outFile" | sed "s/UV/UV_Pool2/g" >> "$outFileTmp"
grep "Pool3.*UV" "$outFile" | sed "s/UV/UV_Pool3/g" >> "$outFileTmp"
#Re-order VIS samples
grep "Pool1.*VIS" "$outFile" | sed "s/VIS/VIS_Pool1/g" >> "$outFileTmp"
grep "Pool2.*VIS" "$outFile" | sed "s/VIS/VIS_Pool2/g" >> "$outFileTmp"
grep "Pool3.*VIS" "$outFile" | sed "s/VIS/VIS_Pool3/g" >> "$outFileTmp"
#Remove excess pool tags
sed -i "s/Pool1_//g" "$outFileTmp"
sed -i "s/Pool2_//g" "$outFileTmp"
sed -i "s/Pool3_//g" "$outFileTmp"
#Separate header for sorting
head -1 "$outFileTmp" > "$outFile"
tail -n +2 "$outFileTmp" | sort -k1 -n -t, >> "$outFile"
#Add method tag to each sample
sed -i "s/$/,$1/" "$outFile"
#Add method tag header
sed -i "s/concordant,$1/concordant,method/" "$outFile"
#Clean up
rm "$outFileTmp"