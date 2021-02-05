#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalSearch_blastp_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub reciprocalSearch_blastp.sh proteomeFasta genomeTranscripts
#Usage Ex: qsub reciprocalSearch_blastp.sh PA42_v4.1_proteins PA42_v4.1_proteins
#Usage Ex: qsub reciprocalSearch_blastp.sh PA42_v4.1_proteins PA42_v4.1_cds
#Usage Ex: qsub reciprocalSearch_blastp.sh PA42_v4.1_proteins PA42_v4.1_transcripts
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98 PA42_v4.1_proteins
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity PA42_v4.1_cds
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity PA42_v4.1_transcripts
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyGenomeTrinity PA42_v4.1_proteins

#Load necessary modules for ND CRC servers
module load bio
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input database for blastp
if [[ "$1" == *assemblyTrinity* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	outputFolder="$inputsPath"/"$1"
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$outputFolder"/decoded_transdecoder/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$outputFolder"/decoded_transdecoder/cdhit.transdecoder.pep
	else 
		inputsPath="$outputFolder"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == *assemblyGenome* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	outputFolder="$inputsPath"/"$1"
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$outputFolder"/decoded_transdecoder/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$outputFolder"/decoded_transdecoder/cdhit.transdecoder.pep
	else 
		inputsPath="$outputFolder"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	outputFolder=$(dirname "$inputsPath")
	inputsPath=$(echo "$outputPath"/decoded_transdecoder/*transdecoder.pep)
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	outputFolder=$(dirname "$inputsPath")
	inputsPath=$(echo "$outputPath"/decoded_transdecoder/*transdecoder.pep)
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Check if DB of transcriptome exsists
inputDB=$(dirname "$inputsPath")
if [ -f "$inputsPath".phr ]; then
	echo "Using exsisting "$inputsPath".phr DB..."
else #Make blastable DB of transcriptome
	#Determine current script location
	currLoc=$(echo $PWD)
	#Move to DB directory
	cd $inputDB
	makeblastdb -in $inputsPath -dbtype prot
	#Move back to script location
	cd $currLoc
fi
#Determine which genome transcript set to use
if [[ "$2" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	reciprocalPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
elif [[ "$2" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	reciprocalPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	reciprocalPath=$(dirname "$reciprocalPath")
	reciprocalPath="$reciprocalPath"/decoded_transdecoder/*.transdecoder.pep
elif [[ "$2" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	reciprocalPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	reciprocalPath=$(dirname "$reciprocalPath")
	reciprocalPath="$reciprocalPath"/decoded_transdecoder/*.transdecoder.pep
else
	#Error message
	echo "Invalid genome fasta entered... exiting!"
	exit 1
fi
#Check if DB of transcriptome exsists
reciprocalDB=$(dirname "$reciprocalPath")
if [ -f "$reciprocalPath".phr ]; then
	echo "Using exsisting "$reciprocalPath".phr DB..."
else #Make blastable DB of transcriptome
	#Move to DB directory
	cd $reciprocalDB
	makeblastdb -in $reciprocalPath -dbtype prot
fi
#Set outputs absolute path
outputFolder="$outputFolder"/reciprocalSearched_blastp_"$2"
#Name output file of inputs
inputOutFile="$outputFolder"/reciprocalSearched_blastp_summary.txt
#Make output directory
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Use blastp to search a database
# and output with outfmt6 header:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Beginning blastp database search..."
blastp -query "$inputsPath" -db "$reciprocalPath" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > "$outputFolder"/blastp.outfmt6
echo "Finished blastp database search!"
#Output run commands to summary file
echo "blastp -query $inputsPath -db $reciprocalPath  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 >" "$outputFolder"/"blastp.outfmt6" >> "$inputOutFile"
#Switch query and search paths for reciprocal search
echo "Beginning reciprocal blastp database search..."
blastp -query "$reciprocalPath" -db "$inputsPath" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > "$outputFolder"/blastp_reciprocal.outfmt6
echo "Finished reciprocal blastp database search!"
#Output run commands to summary file
echo "blastp -query $reciprocalPath -db $inputsPath  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 >" "$outputFolder"/"blastp_reciprocal.outfmt6" >> "$inputOutFile"