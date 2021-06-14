#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateMSA_pep_jobOutput

#Load necessary modules
module load bio

#Retrieve input absolute path
outPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
outPath="$outPath"/"$2"/"$3"
inputsPath="$outPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep
#Retrieve genome reference absolute path
refPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath=$(dirname $refPath)
refPath="$refPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep

#Prepare multiline pep fasta to retrieve seqs
colRefFile="$1"
tmpRef="$colRefFile"_tmpPA42_v4.1.fasta
cat "$refPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRef"

#Loop over all genes in the reference
while IFS= read -r line; do
	#Retrieve selected peptide sequences and convert back to multiline fasta format
	gTag=$(echo $line | sed "s/DASH/-/g" | sed "s/PERIOD/\./g")
	gFile="$outPath"/tmp_pep_allDaphnia_"$line".fasta
	grep "^>$gTag" "$tmpRef" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$gFile"
	
	#Output Status message
	echo "Generating MSA for $line..."
	
	#Prepare multiline pep fasta to retrieve seqs
	tmpSample="$outPath"/tmp_"$line".fasta
	cat "$inputsPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpSample"
	grep "^>$gTag" "$tmpSample" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" >> "$gFile"

	#Create MSA
	mFile="$outPath"/"$line"_pep_allDaphnia_aligned.fasta
	muscle -in "$gFile" -out "$mFile"
	
	#Output status message
	echo "MSA created for $line: $mFile"

	#Clean up
	rm "$outPath"/tmp_"$line".fasta
done < "$colRefFile"

#Clean up
rm "$tmpRef"
rm "$colRefFile"