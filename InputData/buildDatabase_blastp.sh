#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_blastp_jobOutput
#Script to retrieve Uniprot protein sequence database fasta files
#Usage: qsub buildDatabase_blastp.sh dbSelection
#Usage ex: qsub buildDatabase_blastp.sh swissprot
#Usage ex: qsub buildDatabase_blastp.sh ncbi
#Usage ex: qsub buildDatabase_blastp.sh uniprot

#Load necessary modules
module load bio/blast+
if [[ "$1" == "ncbi" ]]; then
	#Retreive ncbi database storage path
	outputPath=$(grep "ncbiDB:" databasePaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$1" == "uniprot" ]]; then
	#Retreive uniprot database storage path
	outputPath=$(grep "uniprotDB:" databasePaths.txt | tr -d " " | sed "s/uniprotDB://g")
elif [[ "$1" == "swissprot" ]]; then
	#Retreive uniprot database storage path
	outputPath=$(grep "swissprotDB:" databasePaths.txt | tr -d " " | sed "s/swissprotDB://g")
else
	#Error message
	echo "Invalid database selection entered (ncbi or uniprot only)... exiting!"
	exit 1
fi
#Make output database directory
outputPath=$(dirname "$outputPath")
mkdir $outputPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	#Error message
	echo "Database directory for $1 files already exsists... exiting!"
	exit 1
fi
#Move to output database directory
cd $outputPath
#Index the database for blastp
dbFileNoEx=$(echo $dbFile | sed 's/\.gz//')
makeblastdb -in $dbFileNoEx -dbtype prot
