#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalSearch_blastp_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub reciprocalSearch_blastp.sh transcriptomeFasta genomeTranscripts
#Usage Ex: qsub reciprocalSearch_blastp.sh PA42_proteins PA42_proteins
#Usage Ex: qsub reciprocalSearch_blastp.sh PA42_proteins PA42_cds
#Usage Ex: qsub reciprocalSearch_blastp.sh PA42_proteins PA42_transcripts
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity PA42_proteins
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity PA42_cds
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity PA42_transcripts

#Load necessary modules for ND CRC servers
module load bio/blast+
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input database for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputDB="$inputsPath"/"$1"
	inputsPath="$inputDB"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
elif [[ "$1" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	inputDB=$(dirname "$inputsPath")
	inputsPath="$inputDB"/PA42.3.0.protein_new.fasta
elif [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputDB=$(dirname "$inputsPath")
	inputsPath="$inputDB"/decoded_transdecoder/PA42.3.0.cds_new.fasta.transdecoder.pep
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputDB=$(dirname "$inputsPath")
	inputsPath="$inputDB"/decoded_transdecoder/PA42.3.0.transcripts_new.fasta.transdecoder.pep
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Check if DB of transcriptome exsists
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
if [[ "$2" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	reciprocalPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	reciprocalDB=$(dirname "$reciprocalPath")
	reciprocalPath="$reciprocalDB"/PA42.3.0.protein_new.fasta
elif [[ "$2" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	reciprocalPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	reciprocalDB=$(dirname "$reciprocalPath")
	reciprocalPath="$reciprocalDB"/decoded_transdecoder/PA42.3.0.cds_new.fasta.transdecoder.pep
elif [[ "$2" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	reciprocalPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	reciprocalDB=$(dirname "$reciprocalPath")
	reciprocalPath="$reciprocalDB"/decoded_transdecoder/PA42.3.0.transcripts_new.fasta.transdecoder.pep
else
	#Error message
	echo "Invalid PA42 fasta entered... exiting!"
	exit 1
fi
#Check if DB of transcriptome exsists
if [ -f "$reciprocalPath".phr ]; then
	echo "Using exsisting "$reciprocalPath".phr DB..."
else #Make blastable DB of transcriptome
	#Move to DB directory
	cd $reciprocalDB
	makeblastdb -in $reciprocalPath -dbtype prot
fi
#Set outputs absolute path
outputFolder="$inputDB"/reciprocalSearched_blastp_"$2"
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