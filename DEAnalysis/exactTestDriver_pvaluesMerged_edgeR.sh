#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash exactTestDriver_pvaluesMerged_edgeR.sh sampleList
#Usage Ex: bash exactTestDriver_pvaluesMerged_edgeR.sh 0.10 Y05 Y023_5 E05 R2 PA Sierra

#Set analysis inputs path
inputsPath="/home/mae/Documents/RNASeq_Workshop_ND"
inFile="$inputsPath"/"geneCounts_cleaned_PA42_v4.1.csv"
#Set FDR cut off
fdrCut=$1

#Set file name for output merged table of pvalues
mergedFile="exactTest_FDR"$fdrCut"_pValues.csv"

#Loop through all input sets of treatments and perform t-test analsysis
for i in "${@:2}"; do
	#Retrieve selected sample column range
	colNumStart=$(($(head -1 "$inFile" | tr "," "\n" | grep -n "$i" | head -1 | cut -d ':' -f1)-1))
	colNumEnd=$(($colNumStart+5))
	#Perform DE analysis using edgeR and output analysis results to a txt file
	Rscript exactTest_edgeR_pvalues.r "$inFile" $colNumStart $colNumEnd $fdrCut $i
done

#Set the initial header for the merged table
echo "gene,$2" > $mergedFile
#Initialize table of merged p-values
tail -n +2 "exactTest"_"$2"_"pvalues.csv" | cut -f 1,4 >> $mergedFile

#Generate merged table of p-values
rm uniq_"$mergedFile"
for i in "${@:3}"; do
	#Output status message
	echo "Merging genotype $i..."
	#Exact test pvalues for the current genotype
	sampleValues="exactTest"_"$i"_"pvalues.csv"
	#Remove headers
	tail -n +2 $sampleValues > tmp_pvalues.txt
	#Update header with current genotype
	sed -i "/^gene,/ s/$/,$i/" $mergedFile
	#Loop over exact test data for the current genotype
	while IFS=, read -r f1 f2 f3 f4
	do
		#Determine if the gene is already in the merged file
		if grep -q "$f1," $mergedFile; then #Matched hits
			sed -i "/^$f1,/ s/$/,$f4/" $mergedFile
		else #Unique hits
			echo "$f1,$f2,$i" >> uniq_"$mergedFile"
		fi
	done < tmp_pvalues.txt
	#Clean up
	rm $sampleValues
done
#Output status message
echo "Merged file of p-values created: $mergedFile"

#Clean up
rm tmp_pvalues.txt