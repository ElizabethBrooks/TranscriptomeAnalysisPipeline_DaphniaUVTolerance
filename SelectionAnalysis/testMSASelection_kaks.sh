#!/bin/bash

#Load necessary module
module load bio

#Set output paths
outPath=$(dirname "$1")
outPath="$outPath"/"$line"MSA_kaks.txt
tOutPath="$outPath"/tmp_"$line"MSA_kaks.txt

#Loop over all genes in the reference
colRefFile="$1"
while IFS= read -r line; do
	#Calculate ka ks values from each MSA
	Rscript testSelection_kaks.r > "$tOutPath"
	#Re format output distance matrix results
	cat "$tOutPath" | sed 's/ \+/,/g' | sed '/^,/d' | sed '/,$/d' > "$outPath"
	#Clean up
	rm "$tOutPath"
done < "$colRefFile"

#Clean up
rm "$colRefFile"