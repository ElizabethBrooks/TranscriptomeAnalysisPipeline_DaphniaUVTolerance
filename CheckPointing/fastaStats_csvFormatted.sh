#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats_csvFormatted.sh fastaFilePaths

#Create header for csv
echo "file, sequences, lines, bytes"

#Initialize counters
counter=1
seqsTotal=0
linesTotal=0
bytesTotal=0

#Loop over input files and retrieve stats
for i in "$@"; do
	seqs=$(grep ">" $i | wc -l)
	lines=$(wc -l $i | cut -d " " -f1)
	bytes=$(wc -c $i | cut -d " " -f1)
	#Output current file name and stats
	echo "file$counter, $seqs, $lines, $bytes"

	#Calculate running total of file stats
	seqsTotal=$(($seqsTotal+$seqs))
	linesTotal=$(($linesTotal+$lines))
	bytesTotal=$(($bytesTotal+$bytes))

	#Increment file tag counter
	counter=$(($counter+1))
done

#Output total file stats
echo "total, $seqsTotal, $linesTotal, $bytesTotal"