#!/bin/bash

#Load necessary module
module load bio

#Set output paths
outPath=$(dirname "$1")
outPath="$outPath"/"$line"MSA_kaks.txt
tOutPath="$outPath"/tmp_"$line"MSA_kaks.txt
fOutPath="$outPath"/tmpFormatted_"$line"MSA_kaks.txt

#Loop over all genes in the reference
colRefFile="$1"
while IFS= read -r line; do
	#Calculate ka ks values from each MSA
	Rscript testSelection_kaks.r > "$tOutPath"
	#Re format output distance matrix results
	cat "$tOutPath" | sed 's/ \+/,/g' | sed '/,$/d' | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE\$/\n\$/g' | grep "^\$ka" | sed 's/NEWLINE/\n/g' | sed 's/\$ka//g' | sed '/^[[:space:]]*$/d' | sed 's/$/,ka/' > "$fOutPath"
	cat "$tOutPath" | sed 's/ \+/,/g' | sed '/,$/d' | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE\$/\n\$/g' | grep "^\$ks" | sed 's/NEWLINE/\n/g' | sed 's/\$ks//g' | sed '/^[[:space:]]*$/d' | sed 's/$/,ks/' >> "$fOutPath"
	cat "$tOutPath" | sed 's/ \+/,/g' | sed '/,$/d' | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE\$/\n\$/g' | grep "^\$vka" | sed 's/NEWLINE/\n/g' | sed 's/\$vka//g' | sed '/^[[:space:]]*$/d' | sed 's/$/,vka/' >> "$fOutPath"
	cat "$tOutPath" | sed 's/ \+/,/g' | sed '/,$/d' | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE\$/\n\$/g' | grep "^\$vks" | sed 's/NEWLINE/\n/g' | sed 's/\$vks//g' | sed '/^[[:space:]]*$/d' | sed 's/$/,vks/' >> "$fOutPath"
	#Re format results with pairwise tags
	echo "sample1,sample2,value,metric" > "$outPath"
	while IFS= read -r line; do
		#echo "Pprocessing line $line"
	  	if $(echo "$line" | grep -q "^,"); then
	    	pastLine="$line"
	  	else
		  	if $(echo "$pastLine" | grep -q "^,"); then
		    	pairTag=$(echo "$pastLine" | sed 's/^,//g' | sed 's/,.*//g' | cut -d ',' -f 1)
		    	echo "$pairTag","$line" >> "$outPath"
		  	fi
	  	fi
	done < "$fOutPath"
	#Clean up
	rm "$tOutPath"
	rm "$fOutPath"
done < "$colRefFile"

#Clean up
rm "$colRefFile"