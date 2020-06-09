#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats_csvFormatted.sh mergedFile fastaFilePaths

#Create header for csv
echo "file,sequences,lines,MB"

#Initialize counters
counter=1
seqsTotal=0
linesTotal=0
mBytesTotal=0

#Loop over input files and retrieve stats
for i; do
	#Calculate current file stats
	seqs=$(grep ">" $i | wc -l)
	lines=$(wc -l $i | cut -d " " -f1)
	mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5)
	#Output file name tags and stats
	echo "file$counter,$seqs,$lines,$mBytes"

	#Calculate running total of un-merged file stats
	if [[ $counter -lt $# ]]; then
		seqsTotal=$(($seqsTotal+$seqs))
		linesTotal=$(($linesTotal+$lines))
		mBytesTotal=$(($mBytesTotal+$mBytes))
	fi
	#Increment counter
	counter=$(($counter+1))
done

#Output total stats
echo "fileTotal,$seqsTotal,$linesTotal,$mBytesTotal"

#Calculate merged file stats
seqsMerged=$(grep ">" $1 | wc -l)
linesMerged=$(wc -l $1 | cut -d " " -f1)
mBytesMerged=$(ls -l --block-size=1MB $1 | cut -d " " -f5)
#Calculate duplicates stats
seqsDup=$(($seqsTotal-$seqsMerged))
linesDup=$(($linesTotal-$linesMerged))
mBytesDup=$(($mBytesTotal-$mBytesMerged))
#Output duplicates stats
echo "fileDuplicates,$seqsDup,$linesDup,$mBytesDup"
