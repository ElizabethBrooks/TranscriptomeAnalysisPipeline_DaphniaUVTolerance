#!/bin/bash
#Script to prepare alignment summary matrices
#Usage: bash alignmentSummary_formatting.sh analysisMethod runNum
#Usage Ex: bash alignmentSummary_formatting.sh tophat2 3 
#Alternate usage Ex: bash alignmentSummary_formatting.sh trimmed_run1E05_genomeTrinity_hisat2_run1 1 

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
#Determine analysis method used
if [[ "$1" == *Trinity* ]]; then
	method="transcriptome"
else
	method="genome"
fi
#Determine alignment software used
if [[ "$1" == *hisat2* ]]; then
	software="hisat2"
else
	software="tophat2"
fi
#Add method, alignment software and run number tag to each sample
sed -i "s/$/,$method,$software,run$2/" "$outFile"
#Add method tag header
sed -i "s/concordant,$method,$software,run$2/concordant,method,software,run/" "$outFile"
#Determine if tophat2 data was input
if [[ "$1" == "tophat2" ]]; then
	#Remove extra columns from tophat2 data
	cut -d, -f2-3 --complement "$outFile" > "$outFileTmp"
	mv "$outFileTmp" "$outFile" 
else
	#Clean up hisat2 tmp file
	rm "$outFileTmp"
fi