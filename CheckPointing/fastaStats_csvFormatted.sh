#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats_csvFormatted.sh fastaFilePaths

#Create header for csv
echo "file, sequences, lines, MB"

#Initialize counters
counter=1
seqsTotal=0
linesTotal=0
mBytesTotal=0

#Loop over input files and retrieve stats
for i in "$@"; do
	#Calculate current file stats
	seqs=$(grep ">" $i | wc -l)
	lines=$(wc -l $i | cut -d " " -f1)
	bytes=$(wc -c $i | cut -d " " -f1)
	mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5)

	#Output file name tags and stats
	if [[ "$i" == *merged* ]]; then
		echo "merged, $seqs, $lines, $mBytes"
	else
		echo "file$counter, $seqs, $lines, $mBytes"
	fi

	#Calculate running total of un-merged file stats
	if [[ $counter != $# ]]; then
		seqsTotal=$(($seqsTotal+$seqs))
		linesTotal=$(($linesTotal+$lines))
		mBytesTotal=$(($mBytesTotal+$mBytes))
	fi

	#Increment file tag counter
	counter=$(($counter+1))
done

#Output total file stats
echo "total, $seqsTotal, $linesTotal, $mBytesTotal"