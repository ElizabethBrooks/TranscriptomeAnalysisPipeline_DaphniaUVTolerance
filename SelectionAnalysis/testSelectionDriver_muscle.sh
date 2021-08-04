#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N testSelection_jobOutput
#Script to generate MSAs for each gene in the reference set of peptide sequences
#Usage: qsub testSelectionDriver_muscle.sh sampleSet variantCallingDir
#Usage ex: qsub testSelectionDriver_muscle.sh sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ

#Load necessary modules
module load bio

#Retrieve input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="$inputsPath"/"$1"/"$2"
inputsPath="$inputsPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep
#Retrieve genome reference absolute path
refPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath=$(dirname $refPath)
refPath="$refPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep

#Set outputs path
outDir=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
outPath="$outDir"/daphniaMSA_PA42_v4.1_pep
#Check if the folder already exists
mkdir "$outPath"
if [ $? -ne 0 ]; then
	echo "The $outPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Set results path
resultsDir="$outDir"/daphniaKaks_PA42_v4.1
#Check if the folder already exists
mkdir "$resultsDir"
if [ $? -ne 0 ]; then
	echo "The $resultsDir directory already exsists... please remove before proceeding."
	exit 1
fi

#Generate a list of all genes in the reference
colRefFile="$outPath"/col1.txt
grep ">" "$refPath" | cut -d" " -f1 | sed "s/-/DASH/g" | sed "s/\./PERIOD/g" | sed "s/>//g" > "$colRefFile"

#Split the reference gene set into segments
#split -l 2000 "$outPath"/col1.txt "$outPath"/colRef

#Generate MSAs for each segment
#for i in "$outPath"/colRef*; do
#	qsub generateMSA_pep_muscle.sh "$i" "$@"
#done

#Prepare reference multiline pep fasta to retrieve seqs
tmpRef="$colRefFile"_tmpPA42_v4.1.fasta
cat "$refPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRef"

#Prepare input multiline pep fasta
tmpSample="$colRefFile"_tmpInput.fasta
cat "$inputsPath" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpSample"

#Save ka ks values to final results file
resultsFile="$resultsDir"/kaksResults.csv
echo "geneID  t  S  N  dNdS  dN  dS" > "$resultsFile"

#Loop over all genes in the reference
outPath=$(dirname $colRefFile)
while IFS= read -r line; do
	#Retrieve selected peptide sequences and convert back to multiline fasta format
	gTag=$(echo $line | sed "s/DASH/-/g" | sed "s/PERIOD/\./g")
	gFile="$outPath"/tmp_pep_allDaphnia_"$line".fasta
	grep "^>$gTag" "$tmpRef" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$gFile"
	
	#Output Status message
	echo "Generating MSA for $gTag..."

	#Prepare multiline pep fasta to retrieve seqs
	grep "^>$gTag" "$tmpSample" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" >> "$gFile"

	#Create MSA
	mFile="$outPath"/"$gTag"_pep_allDaphnia_aligned.fasta
	muscle -in "$gFile" -out "$mFile"
	
	#Output status message
	echo "MSA created for $gTag: $mFile"

	#Clean up
	rm "$gFile"

	#Generate and save ka ks values
	bash testSelection_pal2nalCodeml.sh "$gTag" "$1" "$2"
done < "$colRefFile"

#Fix formatting of the results file
finalResults="$resultsDir"/PA42_v4.1_Olympics_kaksResults.csv
cat "$resultsFile" | sed "s/  /,/g" | sed "s/=,/=/g" | sed "s/ //g" | sed "s/pairwisecomparison,codonfrequencies\:F3x4\./NA,NA,NA,NA,NA,NA/g" | sed "s/dN\/dS=//g" | sed "s/dN=//g" | sed "s/dS=//g" | sed "s/t=//g" | sed "s/S=//g" | sed "s/N=//g" > "$finalResults"

#Clean up
rm "$tmpSample"
rm "$tmpRef"
rm "$colRefFile"
rm -r "$outPath"
rm "$resultsFile"