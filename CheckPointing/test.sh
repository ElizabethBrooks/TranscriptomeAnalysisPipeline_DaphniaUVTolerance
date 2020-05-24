#!/bin/bash
#Script to merge multiple fasta files and retain only
#the specified unique data (by sequence, name, or sequenceAndName)
#Usage: bash test.sh mergeBy fastaFilePaths
#Usage Ex: bash test.sh sequence ~/trimmed_run1/Trinity.fasta ~/trimmed_run2/Trinity.fasta

#Check for input arguments of fasta files
if [ $# -eq 0 ]; then
   	echo "No fasta files input... exiting!"
   	exit 1
fi

#Store list of input fasta file paths after
# skipping the first input argument of mergeBy
for i in "${@:2}"; do
    fastaList="$fastaList$i "
done

#Set output file paths
outputFastaFile="./Trinity.fasta"
summaryFile="./mergedFasta_summary.txt"



#Write fasta stats to the summary file
bash fastaStats.sh $fastaList $outputFastaFile > $summaryFile

#Write fasta stats to the csv formatted summary file
summaryFileCSV=$(echo "$summaryFile" | sed 's/\.txt/\.csv/g')
bash fastaStats_csvFormatted.sh $fastaList $outputFastaFile > $summaryFileCSV

#Plot fasta stats from summary file
Rscript fastaStats_barPlot.r $summaryFileCSV $1

#Clean up
rm Rplots.pdf