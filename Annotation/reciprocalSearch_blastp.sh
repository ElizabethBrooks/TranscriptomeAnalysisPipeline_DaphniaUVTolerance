#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalSearch_blastp_jobOutput
#$ -pe smp 8

# Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage 1: qsub reciprocalSearch_blastp.sh searchTag proteomeFastaQuery
#Usage 2: qsub reciprocalSearch_blastp.sh searchTag proteomeFastaDB proteomeFastaQuery
#Usage ex: qsub reciprocalSearch_blastp.sh Dpulex_Dmelanogaster /afs/crc.nd.edu/group/pfrenderlab/mendel/OutgroupGenomes/ncbi_dataset_Drosophila_melanogaster/ncbi_dataset/data/GCF_000001215.4/protein.faa

# Load necessary modules for ND CRC servers
module load bio/2.0

# Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
elif [ $# -eq 1 ]; then # use PA42 protein seqs
	# Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
	#Retrieve database absolute path for querying
	reciprocalPath="$1"
elif [ $# -eq 2 ]; then # expect two paths to protein seqs
	# Retrieve genome reference absolute path for querying
	inputsPath="$1"
	#Retrieve database absolute path for querying
	reciprocalPath="$2"
fi

# set outputs path
outputFolder=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")

# set outputs absolute folder name
searchTag="$1"
outputFolder=$outputFolder"/reciprocalSearched_blastp_"$searchTag
#Name output file of inputs
inputOutFile=$outputFolder"/reciprocalSearched_blastp_summary.txt"

#Make output directory
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"

# Check if first DB of seqs exsists
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

# Check if second DB of seqs exsists
reciprocalDB=$(dirname "$reciprocalPath")
if [ -f "$reciprocalPath".phr ]; then
	echo "Using exsisting "$reciprocalPath".phr DB..."
else #Make blastable DB of seqs
	#Move to DB directory
	cd $reciprocalDB
	makeblastdb -in $reciprocalPath -dbtype prot
fi

#Use blastp to search a database
# and output with outfmt6 header:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Beginning blastp database search..."
blastp -query "$inputsPath" -db "$reciprocalPath" -outfmt 6 -evalue 0.01 -num_threads 8 > "$outputFolder"/blastp.outfmt6
echo "Finished blastp database search!"
#Output run commands to summary file
echo "blastp -query $inputsPath -db $reciprocalPath  -outfmt 6 -evalue 0.01 -num_threads 8 >" "$outputFolder"/"blastp.outfmt6" >> "$inputOutFile"
#Switch query and search paths for reciprocal search
echo "Beginning reciprocal blastp database search..."
blastp -query "$reciprocalPath" -db "$inputsPath" -outfmt 6 -evalue 0.01 -num_threads 8 > "$outputFolder"/blastp_reciprocal.outfmt6
echo "Finished reciprocal blastp database search!"
#Output run commands to summary file
echo "blastp -query $reciprocalPath -db $inputsPath -outfmt 6 -evalue 0.01 -num_threads 8 >" "$outputFolder"/"blastp_reciprocal.outfmt6" >> "$inputOutFile"