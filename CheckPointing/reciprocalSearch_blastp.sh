#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalSearch_blastp_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub reciprocalSearch_blastp.sh transcriptomeFasta
#Usage Ex: qsub reciprocalSearch_blastp.sh trimmed_run1E05_assemblyTrinity

#Load necessary modules for ND CRC servers
module load bio/blast+
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve genome reference absolute path for querying
dbPath=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
#Determine input database for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/reciprocalSearchedTranscripts_blastp
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
inputOutFile="$outputFolder"/reciprocalSearchedTranscripts_blastp_"$1"_summary.txt
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