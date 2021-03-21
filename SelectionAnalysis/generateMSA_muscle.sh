#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateMSA_jobOutput

#Load necessary modules
module load bio

#Retrieve inputs path
inputPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
refPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
#Set outputs path
outPath=$(dirname "$1")

#Loop over all genes in the reference
colRefFile="$1"
while IFS= read -r line; do
	#Prepare multiline cds fasta to retrieve seqs
	cat "$refPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$outPath"/tmpPA42_v4.1.fasta
	#Retrieve selected coding sequences and convert back to multiline fasta format
	gTag=$line"-CDS"
	gFile="$outPath"/"$line"_cds_allDaphnia.fasta
	grep "^>$gTag" "$outPath"/tmpPA42_v4.1.fasta | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$gFile"
	
	#Output Status message
	echo "Generating MSA for $line..."
	#Loop over all input samples
	for i in "${@:2}"; do
		#Prepare multiline cds fasta to retrieve seqs
		cat "$inputPath"/sortedCoordinate_samtoolsHisat2_run1"$i"_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$outPath"/tmp"$i".fasta
		#Retrieve gene associated with RBHB
		pTag=$line"-pep"
		rbhTag=$(grep "$pTag" "$inputPath"/sortedCoordinate_samtoolsHisat2_run1"$i"_assemblyPA42_v4.1Trinity/reciprocalSearched_blastp_PA42_v4.1_proteins/blastp_RBH.txt | cut -d ',' -f 1)
		#Retrieve selected coding sequences and convert back to multiline fasta format
		sTag="$i"_"$rbhTag"
		grep "^>$rbhTag" "$outPath"/tmp"$i".fasta | sed 's/NEWLINE/\n/g' | sed "s/^>$rbhTag.*/>$sTag/g" >> "$gFile"
	done

	#Create MSA
	mFile="$outPath"/"$line"_cds_allDaphnia_aligned.fasta
	muscle -in "$gFile" -out "$mFile"
	#Output status message
	echo "$line MSA created: $mFile"

	#Clean up
	rm "$outPath"/tmp*.fasta
	rm "$gFile"
done < "$colRefFile"

#Clean up
rm "$colRefFile"