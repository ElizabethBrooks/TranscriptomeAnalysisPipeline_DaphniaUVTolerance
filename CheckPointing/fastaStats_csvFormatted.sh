#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats_csvFormatted.sh fastaFilePaths

#Create header for csv
echo "file, sequences, lines, MB"

#Initialize counters
counter=1

#Loop over input files and retrieve stats
for i in "$@"; do
	echo "Calculating stats for file $i..."
	#Calculate current file stats
	seqs=$(grep ">" $i | wc -l)
	lines=$(wc -l $i | cut -d " " -f1)
	bytes=$(wc -c $i | cut -d " " -f1)
	mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5)

	#Output file name tags and stats
	echo "file$counter, $seqs, $lines, $mBytes"

	#Increment file tag counter
	counter=$(($counter+1))
done
