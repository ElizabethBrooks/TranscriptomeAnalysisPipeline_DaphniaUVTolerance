#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalSearch_blastp_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub reciprocalSearch_blastp.sh transcriptomeFasta genomeTranscripts
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
#Determine current script location
currLoc=$(echo $PWD)
#Determine input database for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/reciprocalSearched_blastp_"$2"
	#Make blastable DB of transcriptome
	cd $inputsPath
	inputDB=Trinity.fasta.transdecoder.pep
	makeblastdb -in $inputDB -dbtype prot
	inputsPath="$inputsPath"/"$inputDB"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Move back to script location
cd $currLoc
#Determine which genome transcript set to use
if [[ "$2" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	dbPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	#Make blastable DB of transcripts
	cd $dbPath
	makeblastdb -in $dbPath -dbtype prot
elif [[ "$2" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	dbPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputDB=$(dirname "$dbPath")
	dbPath="$inputDB"/decoded_transdecoder/PA42.3.0.cds_new.fasta.transdecoder.pep
	#Make blastable DB of transcripts
	cd $dbPath
	makeblastdb -in $dbPath -dbtype prot
elif [[ "$2" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	dbPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputDB=$(dirname "$dbPath")
	dbPath="$inputDB"/decoded_transdecoder/PA42.3.0.transcripts_new.fasta.transdecoder.pep
	#Make blastable DB of transcripts
	cd $dbPath
	makeblastdb -in $dbPath -dbtype prot
else
	#Error message
	echo "Invalid PA42 fasta entered... exiting!"
	exit 1
fi
#Make output directory
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/reciprocalSearched_blastp_"$1"_summary.txt
#Use blastp to search a database
# and output with outfmt6 header:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Beginning blastp database search..."
blastp -query "$inputsPath" -db "$dbPath" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp.outfmt6
echo "Finished blastp database search!"
#Output run commands to summary file
echo "blastp -query $inputsPath -db $dbPath  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp.outfmt6" >> "$inputOutFile"
#Switch query and search paths for reciprocal search
echo "Beginning reciprocal blastp database search..."
blastp -query "$dbPath" -db "$inputsPath" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp_reciprocal.outfmt6
echo "Finished reciprocal blastp database search!"
#Output run commands to summary file
echo "blastp -query $dbPath -db $inputsPath  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp_reciprocal.outfmt6" >> "$inputOutFile"