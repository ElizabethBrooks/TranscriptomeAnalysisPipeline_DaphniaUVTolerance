#!/bin/bash
#Script to retrieve stats for an input list of files
#Usage: bash fastaStats.sh fileList

#Check number of sequences
echo "Number of sequences:"
for i in "$@"; do
	seqs=$(grep ">" $i | wc -l); echo "$seqs $i"
done
#Check number of lines
echo "Number of lines:"
for i in "$@"; do
	wc -l $i
done
#Check file sizes
echo "File sizes (bytes):"
for i in "$@"; do
	wc -c $i
done