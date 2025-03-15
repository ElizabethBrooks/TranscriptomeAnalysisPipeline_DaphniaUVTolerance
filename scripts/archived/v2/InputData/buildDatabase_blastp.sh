#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N buildDatabase_blastp_jobOutput
#Script to retrieve Uniprot protein sequence database fasta files
#Usage: qsub buildDatabase_blastp.sh dbSelection
#Usage ex: qsub buildDatabase_blastp.sh swissprot
#Usage ex: qsub buildDatabase_blastp.sh ncbi
#Usage ex: qsub buildDatabase_blastp.sh uniprot

#Load necessary modules
module load bio/blast+
#Determine database selection
if [[ "$1" == "ncbi" ]]; then
	#Retreive ncbi database storage path
	dbFile=$(grep "ncbiDB:" databasePaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$1" == "uniprot" ]]; then
	#Retreive uniprot database storage path
	dbFile=$(grep "uniprotDB:" databasePaths.txt | tr -d " " | sed "s/uniprotDB://g")
elif [[ "$1" == "swissprot" ]]; then
	#Retreive uniprot database storage path
	dbFile=$(grep "trinotateSwissprotDB:" databasePaths.txt | tr -d " " | sed "s/trinotateSwissprotDB://g")
else
	#Error message
	echo "Invalid database selection entered (ncbi or uniprot only)... exiting!"
	exit 1
fi
#Move to output database directory
outputPath=$(dirname "$dbFile")
cd $outputPath
#Index the database for blastp
makeblastdb -in $dbFile -dbtype prot
