#!/bin/bash
#Script to prepare alignment summary matrices
#Usage: bash alignmentSummary_formatting.sh alignedFolder analysisMethod runNum
#Usage Ex: bash alignmentSummary_formatting.sh aligned_tophat2_run3 tophat2 run3

#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory and file
outDir="$outputsPath"/AlignmentsAnalyzed
outFile="$outDir"/alignmentSummarized_"$2"_"$3"_formatted.csv
outFileTmp="$outDir"/alignmentSummarized_"$2"_"$3"_tmp.csv
#Remove percent symbols
sed "s/%//g" "$1" > "$outFile"
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
grep "Pool1.*UV" "$1" | sed "s/UV/UV_Pool1/g" >> "$outFileTmp"
grep "Pool2.*UV" "$1" | sed "s/UV/UV_Pool2/g" >> "$outFileTmp"
grep "Pool3.*UV" "$1" | sed "s/UV/UV_Pool3/g" >> "$outFileTmp"
#Re-order VIS samples
grep "Pool1.*VIS" "$1" | sed "s/VIS/VIS_Pool1/g" >> "$outFileTmp"
grep "Pool2.*VIS" "$1" | sed "s/VIS/VIS_Pool2/g" >> "$outFileTmp"
grep "Pool3.*VIS" "$1" | sed "s/VIS/VIS_Pool3/g" >> "$outFileTmp"
#Remove excess pool tags
sed -i "s/Pool1_//g" "$outFileTmp"
sed -i "s/Pool2_//g" "$outFileTmp"
sed -i "s/Pool3_//g" "$outFileTmp"
#Separate header for sorting
head -1 "$outFileTmp" > "$outFile"
tail -n +2 "$outFileTmp" | sort -k1 -n -t, >> "$outFile"
#Add method tag to each sample
sed -i 's/$/,"$2"/' "$outFile"
#Add method tag header
sed -i 's/concordant,"$2"/concordant,method/' "$outFile"
#Clean up
rm "$outFileTmp"