#!/bin/bash
#Script to generate MSAs for each gene in the reference set of peptide sequences
#Usage: bash generateMSADriver_pep_muscle.sh sampleSet
#Usage ex: bash generateMSADriver_pep_muscle.sh sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ

#Retrieve input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="$inputsPath"/"$1"/"$2"
inputsPath="$inputsPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep
#Retrieve genome reference absolute path
refPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath=$(dirname $refPath)
refPath="$refPath"/"decoded_transdecoder"/transcripts_cufflinks.fa.transdecoder.pep

#Set outputs path
outPath=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
outPath="$outPath"/daphniaMSA_PA42_v4.1_pep

#Check if the folder already exists
mkdir "$outPath"
if [ $? -ne 0 ]; then
	echo "The $outPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Generate a list of all genes in the reference
grep ">" "$refPath" | cut -d" " -f1 | sed "s/-/DASH/g" | sed "s/\./PERIOD/g" | sed "s/>//g" > "$outPath"/col1.txt

#Split the reference gene set into segments
split -l 2000 "$outPath"/col1.txt "$outPath"/colRef

#Generate MSAs for each segment
for i in "$outPath"/colRef*; do
	qsub generateMSA_pep_muscle.sh "$i" "$@"
done

#Clean up
rm "$outPath"/col1.txt