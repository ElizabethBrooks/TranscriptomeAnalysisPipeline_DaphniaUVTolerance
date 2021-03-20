#!/bin/bash

#Retrieve inputs path
inputPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
refPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")

#Loop over all genes in the reference
grep ">" "$refPath" | sed "s/-CDS//g" | sed "s/>//g" > col1.txt
while IFS= read -r line; do
	#Prepare multiline cds fasta to retrieve seqs
	cat "$refPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA42_v4.1.fasta
	#Retrieve selected coding sequences and convert back to multiline fasta format
	gTag=$line"-CDS"
	gFile=$line"_cds_allDaphnia.fasta"
	grep "^>$gTag" tmpPA42_v4.1.fasta | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$gFile"
	#Loop over all input samples
	for i in "$@"; do
		#Prepare multiline cds fasta to retrieve seqs
		cat "$inputPath"/sortedCoordinate_samtoolsHisat2_run1"$i"_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmp"$i".fasta

		#Retrieve gene associated with RBHB
		pTag=$line"-pep"
		rbhTag=$(grep "$pTag" "$inputPath"/sortedCoordinate_samtoolsHisat2_run1"$i"_assemblyPA42_v4.1Trinity/reciprocalSearched_blastp_PA42_v4.1_proteins/blastp_RBH.txt | cut -d ',' -f 1)

		#Retrieve selected coding sequences and convert back to multiline fasta format
		grep "^>$rbhTag" tmp"$i".fasta | sed 's/NEWLINE/\n/g' | sed "s/^>$rbhTag.*/>$i_$rbhTag/g" >> "$gFile"

		#Create MSA
		mFile=$line"_cds_allDaphnia_aligned.fasta"
		muscle -in "$gFile" -out "$mFile"

		#Clean up
		rm tmp*.fasta
		rm "$gFile"
	done
done < col1.txt
#Clean up
rm col1.txt