#!/bin/bash
#Script to merge multiple fasta files and retain only
#the specified unique data (by sequence, name, or sequenceAndName)
#Usage: bash test.sh mergeBy fastaFilePaths
#Usage Ex: bash test.sh sequence ~/trimmed_run1/Trinity.fasta ~/trimmed_run2/Trinity.fasta

module load R

#Check for input arguments of fasta files
if [ $# -lt 2 ]; then
   	echo "Not enough fasta files input... exiting!"
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

#Merge the set of fasta files
bash fastaMerge_driver.sh $1 $outputFastaFile $fastaList

#Write fasta stats to the summary file
bash fastaStats.sh $outputFastaFile $fastaList > $summaryFile

#Write fasta stats to the csv formatted summary file
summaryFileCSV=$(echo "$summaryFile" | sed 's/\.txt/\.csv/g')
bash fastaStats_csvFormatted.sh $outputFastaFile $fastaList > $summaryFileCSV

#Plot fasta stats from summary file
Rscript fastaStats_barPlot.r $1 $summaryFileCSV

#Clean up
rm Rplots.pdf