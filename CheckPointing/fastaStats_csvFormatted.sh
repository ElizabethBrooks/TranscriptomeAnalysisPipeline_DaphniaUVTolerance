#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats_csvFormatted.sh fastaFilePaths

#Create header for csv
echo "file, sequences, lines, bytes"

#Initialize file tag counter
counter=1

#Loop over input files and retrieve stats
for i in "$@"; do
	seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
	lines=$(wc -l $i)
	bytes=$(wc -c $i)
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
echo "totalStats, $seqsTotal, $linesTotal, $bytesTotal"