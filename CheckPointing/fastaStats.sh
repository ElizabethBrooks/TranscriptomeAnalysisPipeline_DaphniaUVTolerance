#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats.sh fileList

#Initialize counters
seqsTotal=0
linesTotal=0
mBytesTotal=0

#Check number of sequences
printf "\nNumber of sequences:"
for i in "$@"; do
	seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
	#Calculate running total of un-merged file stats
	if [[ "$i" != *merged* ]]; then
		seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
		seqsTotal=$(($seqsTotal+$seqs))
	fi
done
echo "$seqsTotal un-merged file total"

#Check number of lines
printf "\nNumber of lines:"
for i in "$@"; do
	wc -l $i
	#Calculate running total of un-merged file stats
	if [[ "$i" != *merged* ]]; then
		lines=$(wc -l $i)
		linesTotal=$(($linesTotal+$lines))
	fi
done
echo "$linesTotal un-merged file total"

#Check file sizes
printf "\nFile sizes (bytes):"
for i in "$@"; do
	mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5); echo "$mBytes $i"
	#Calculate running total of un-merged file stats
	if [[ "$i" != *merged* ]]; then
		mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5)
		mBytesTotal=$(($mBytesTotal+$mBytes))
	fi
done
echo "$mBytesTotal un-merged file total"
