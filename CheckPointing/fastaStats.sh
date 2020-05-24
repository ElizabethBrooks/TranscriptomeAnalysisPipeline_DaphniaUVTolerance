#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats.sh fileList

#Initialize counters
counter=1
seqsTotal=0
linesTotal=0
mBytesTotal=0
lines=0
bytes=0

#Check number of sequences
printf "\nNumber of sequences:"
for i in "$@"; do
	if [[ $counter -lt $# ]]; then
		#Calculate running total of un-merged file stats
		seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
		seqsTotal=$(($seqsTotal+$seqs))
	else
		#Calculate current file stats
		seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
	fi
	#Increment counter
	counter=$(($counter+1))
done
echo "$seqsTotal un-merged file total"

#Re-set counter
counter=1
#Check number of lines
printf "\nNumber of lines:"
for i in "$@"; do
	if [[ $counter -lt $# ]]; then
		#Calculate running total of un-merged file stats
		lines=$(wc -l $i)
		linesTotal=$(($linesTotal+$lines))
	else
		#Calculate current file stats
		wc -l $i
	fi
	#Increment counter
	counter=$(($counter+1))
done
echo "$linesTotal un-merged file total"

#Re-set counter
counter=1
#Check file sizes
printf "\nFile sizes (bytes):"
for i in "$@"; do
	if [[ $counter -lt $# ]]; then
		#Calculate running total of un-merged file stats
		mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5)
		mBytesTotal=$(($mBytesTotal+$mBytes))
	else
		#Calculate current file stats
		mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5); echo "$mBytes $i"
	fi
	#Increment counter
	counter=$(($counter+1))
done
echo "$mBytesTotal un-merged file total"
