#!/bin/bash
#Script to remove duplicate reads from sorted bam files

#Loop over MapQ filtered bam files
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/filteredMapQ.bam; do 
	echo "Processing file $f"
	path=$(dirname $f)
	file=$(basename $f | sed "s/.bam//g")
	#Remove duplicate reads
	samtools markdup -r "$f" "$path"/"$file"_noDups.bam
done