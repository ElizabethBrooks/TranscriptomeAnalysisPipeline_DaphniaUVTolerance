#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateMSA_pep_jobOutput
#Usage: qsub generateMSA_pep_muscle.sh sequenceSubset sampleSet variantCallingDir
#Usage ex: qsub generateMSA_pep_muscle.sh A sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ

#Load necessary modules
module load bio

#Retrieve input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="$inputsPath"/"$2"/"$3"
inputsPath="$inputsPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep
#Retrieve genome reference absolute path
refPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath=$(dirname $refPath)
refPath="$refPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep

#Prepare reference multiline pep fasta to retrieve seqs
colRefFile="$1"
tmpRef="$colRefFile"_tmpPA42_v4.1.fasta
cat "$refPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRef"

#Prepare input multiline pep fasta
tmpSample="$outPath"/"$colRefFile"_tmpInput.fasta
cat "$inputsPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpSample"

#Loop over all genes in the reference
outPath=$(dirname $colRefFile)
while IFS= read -r line; do
	#Retrieve selected peptide sequences and convert back to multiline fasta format
	gTag=$(echo $line | sed "s/DASH/-/g" | sed "s/PERIOD/\./g")
	gFile="$outPath"/tmp_pep_allDaphnia_"$line".fasta
	grep "^>$gTag" "$tmpRef" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$gFile"
	
	#Output Status message
	echo "Generating MSA for $line..."

	#Prepare multiline pep fasta to retrieve seqs
	grep "^>$gTag" "$tmpSample" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" >> "$gFile"

	#Create MSA
	mFile="$outPath"/"$gTag"_pep_allDaphnia_aligned.fasta
	muscle -in "$gFile" -out "$mFile"
	
	#Output status message
	echo "MSA created for $line: $mFile"

	#Generate ka ks values
	bash testSelection_pal2nalCodeml.sh "$gTag" "$2" "$3"
done < "$colRefFile"

#Clean up
rm "$tmpSample"
rm "$tmpRef"
rm "$colRefFile"