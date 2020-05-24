#!/bin/bash
#Script to merge multiple fasta files and retain only
#the specified unique data (by sequence, name, or sequenceAndName)
#Usage: bash mergeFasta.sh mergeBy fastaFilePaths
#Usage Ex: bash mergeFasta.sh sequence ~/trimmed_run1/Trinity.fasta ~/trimmed_run2/Trinity.fasta

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

#Merge a set of fasta files
echo "Beginning fasta file merging..."
#Determine which method to merge fasta files by
if [[ "$1" == sequence ]]; then
	#Sequence identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		!(seen[seq]++){ print ">" $0 }' $fastaList > $outputFastaFile
elif [[ "$1" == name ]]; then
	#First part of sequence name identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		{ key=substr(name,1,index(s,"|")) }
		!(seen[key]++){ print ">" $0 }' $fastaList > $outputFastaFile
elif [[ "$1" == sequenceAndName ]]; then
	#Sequence name and sequence identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		!(seen[name,seq]++){ print ">" $0 }' $fastaList > $outputFastaFile
else
	echo "Selected merge format for fasta files not valid... exiting!"
	exit 1
fi
echo "Fasta file merging complete!"

#Write fasta stats to the summary file
bash fastaStats.sh $fastaList $outputFastaFile > $summaryFile