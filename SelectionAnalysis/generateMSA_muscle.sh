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

#Prepare multiline cds fasta to retrieve seqs
tmpRef="$colRefFile"_tmpPA42_v4.1.fasta
cat "$refPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRef"

#Loop over all genes in the reference
colRefFile="$1"
while IFS= read -r line; do
	#Retrieve selected coding sequences and convert back to multiline fasta format
	gTag=$line"-CDS"
	gFile="$outPath"/tmp_cds_allDaphnia_"$line".fasta
	grep "^>$gTag" "$tmpRef" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$gFile"
	
	#Output Status message
	echo "Generating MSA for $line..."
	#Loop over all input samples
	for i in "${@:2}"; do
		#Prepare multiline cds fasta to retrieve seqs
		tmpSample="$outPath"/tmp_"$i"_"$line".fasta
		inSampleC="$inputPath"/sortedCoordinate_samtoolsHisat2_run1"$i"_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds
		cat "$inSampleC" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpSample"
		#Retrieve gene associated with RBHB
		pTag=$line"-pep"
		inSampleP="$inputPath"/sortedCoordinate_samtoolsHisat2_run1"$i"_assemblyPA42_v4.1Trinity/reciprocalSearched_blastp_PA42_v4.1_proteins/blastp_RBH.txt
		rbhTag=$(grep "$pTag" "$inSampleP" | cut -d ',' -f 1)
		#Retrieve selected coding sequences and convert back to multiline fasta format
		sTag="$i"_"$rbhTag"
		grep "^>$rbhTag" "$tmpSample" | sed 's/NEWLINE/\n/g' | sed "s/^>$rbhTag.*/>$sTag/g" >> "$gFile"
	done

	#Create MSA
	mFile="$outPath"/"$line"_cds_allDaphnia_aligned.fasta
	muscle -in "$gFile" -out "$mFile"
	#Output status message
	echo "MSA created for $line: $mFile"

	#Clean up
	rm "$outPath"/tmp_*_"$line".fasta
done < "$colRefFile"

#Clean up
rm "$tmpRef"
rm "$colRefFile"