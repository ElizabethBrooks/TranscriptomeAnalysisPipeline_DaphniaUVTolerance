#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_blastp_jobOutput
#Script to retrieve Uniprot protein sequence database fasta files
#Usage: qsub buildDatabase_blastp.sh dbSelection
#Usage ex: qsub buildDatabase_blastp.sh ncbi
#Note that the NCBI NR DB may be downloaded with wget
#ex: wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#Note that uniprot databases may be downloaded from the UniprotKB search page, or wget
#ex: wget https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:arthropod
#or ex: wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz

#Load necessary modules
module load bio/blast+
if [[ "$1" == "ncbi" ]]; then
	#Set slected database to ncbi
	dbAddress="ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
	#Retreive ncbi database storage path
	outputPath=$(grep "ncbiDB:" inputPaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$1" == "uniprot" ]]; then
	#Set slected database to uniprot
	dbAddress="ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
	#Retreive uniprot database storage path
	outputPath=$(grep "uniprotDB:" inputPaths.txt | tr -d " " | sed "s/uniprotDB://g")
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
#Retrieve selected input database
wget $dbAddress
#Extract the database
dbFile=$(basename $dbAddress)
gunzip -v $dbFile
#Index the database for blastp
dbFileNoEx=$(echo $dbFile | sed 's/\.gz//')
makeblastdb -in $dbFileNoEx -dbtype prot
