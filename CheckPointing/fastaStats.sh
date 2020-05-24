#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats.sh fileList

#Initialize counters
counter=1
seqsTotal=0
linesTotal=0
mBytesTotal=0

#Output file list with tags
printf "File list:\n"
for i in "$@"; do
	echo "file$counter $i"
	#Increment counter
	counter=$(($counter+1))
done
#Re-set counter
counter=1

#Check number of sequences
printf "\nNumber of sequences:\n"
for i in "$@"; do
	#Calculate current file stats
	seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
	#Calculate running total of un-merged file stats
	if [[ $counter -gt 1 ]]; then
		seqsTotal=$(($seqsTotal+$seqs))
	fi
	#Increment counter
	counter=$(($counter+1))
done
echo "$seqsTotal un-merged file total"
#Re-set counter
counter=1

#Check number of lines
printf "\nNumber of lines:\n"
for i in "$@"; do
	#Calculate current file stats
	lines=$(wc -l $i | awk '{print $1}'); echo "$lines $i"
	#Calculate running total of un-merged file stats
	if [[ $counter -gt 1 ]]; then
		linesTotal=$(($linesTotal+$lines))
	fi
	#Increment counter
	counter=$(($counter+1))
done
echo "$linesTotal un-merged file total"
#Re-set counter
counter=1

#Check file sizes
printf "\nFile sizes (MB):\n"
for i in "$@"; do
	#Calculate current file stats
	mBytes=$(ls -l --block-size=1MB $i | cut -d " " -f5); echo "$mBytes $i"
	#Calculate running total of un-merged file stats
	if [[ $counter -gt 1 ]]; then
		mBytesTotal=$(($mBytesTotal+$mBytes))
	fi
	#Increment counter
	counter=$(($counter+1))
done
echo "$mBytesTotal un-merged file total"
