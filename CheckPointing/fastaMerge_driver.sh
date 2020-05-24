#!/bin/bash
#Script to perform merge mergedfasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: bash fastaMerge_driver.sh mergeBy mergeFileName sortedFolderList

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder inputs supplied... exiting"
   	exit 1
fi

#Store list of input fasta file paths after
# skipping the first two input arguments
for i in "${@:3}"; do
    fastaList="$fastaList$i "
done

#Merge the set of fasta files
echo "Beginning fasta file merging..."
#Determine which method to merge fasta files by
if [[ "$1" == sequence ]]; then
	#Sequence identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		!(seen[seq]++){ print ">" $0 }' $fastaList > $2
elif [[ "$1" == name ]]; then
	#First part of sequence name identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		{ key=substr(name,1,index(s,"|")) }
		!(seen[key]++){ print ">" $0 }' $fastaList > $2
elif [[ "$1" == sequenceAndName ]]; then
	#Sequence name and sequence identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		!(seen[name,seq]++){ print ">" $0 }' $fastaList > $2
else
	echo "Selected merge format for fasta files not valid... exiting!"
	exit 1
fi
echo "Fasta file merged to $2"
