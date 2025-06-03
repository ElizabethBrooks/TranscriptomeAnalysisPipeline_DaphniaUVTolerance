#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash exactTestDriver_pvalues_edgeR.sh sampleList
#Usage Ex: bash exactTestDriver_pvalues_edgeR.sh 0.10 Y05 Y023_5 E05 R2 PA Sierra

#Set analysis inputs path
inputsPath="/home/mae/Documents/RNASeq_Workshop_ND"
inFile="$inputsPath"/"geneCounts_cleaned_PA42_v4.1.csv"

#Loop through all input sets of treatments and perform t-test analsysis
for i in "${@:2}"; do
	#Retrieve selected sample column range
	colNumStart=$(($(head -1 "$inFile" | tr "," "\n" | grep -n "$i" | head -1 | cut -d ':' -f1)-1))
	colNumEnd=$(($colNumStart+5))
	#Perform DE analysis using edgeR and output analysis results to a txt file
	Rscript exactTest_edgeR_pvalues.r "$inFile" $colNumStart $colNumEnd $i
done
