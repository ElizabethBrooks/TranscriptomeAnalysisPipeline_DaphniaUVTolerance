#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N getDatabase_jobOutput
#Script to retrieve Uniprot protein sequence database fasta files
#Usage: qsub getDatabase.sh dbSelection
#Usage ex: qsub getDatabase.sh swissprot
#Usage ex: qsub getDatabase.sh ncbi
#Note that the NCBI NR DB may be downloaded with wget
#ex: wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#Note that uniprot databases may be downloaded from the UniprotKB search page, or wget
#ex: wget https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:arthropod
#or ex: wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz

#Determine which database was selected
if [[ "$1" == "ncbi" ]]; then
	#Set slected database to ncbi
	dbAddress="ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
	#Retreive ncbi database storage path
	outputPath=$(grep "ncbiDB:" databasePaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$1" == "uniprot" ]]; then
	#Set slected database to uniprot
	dbAddress="ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
	#Retreive uniprot database storage path
	outputPath=$(grep "uniprotDB:" databasePaths.txt | tr -d " " | sed "s/uniprotDB://g")
elif [[ "$1" == "swissprot" ]]; then
	#Retrieve Trinotate software path
	softsPath=$(grep "trinotatePackage:" softwarePaths.txt | tr -d " " | sed "s/trinotatePackage://g")
	#Retreive uniprot database storage path
	outputPath=$(grep "swissprotDB:" databasePaths.txt | tr -d " " | sed "s/swissprotDB://g")
elif [[ "$1" == "panther" ]]; then
	#Set slected database to panther
	dbAddress="ftp://ftp.pantherdb.org/panther_library/14.1/PANTHER14.1_hmmscoring.tgz"
	#Retreive uniprot database storage path
	outputPath=$(grep "pantherDB:" databasePaths.txt | tr -d " " | sed "s/pantherDB://g")	
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
if [[ "$1" == "ncbi" || "$1" == "uniprot" ]]
	#Retrieve selected input database
	wget $dbAddress
	#Extract the database
	dbFile=$(basename $dbAddress)
	gunzip -v $dbFile
elif [[ "$1" == "panther" ]]; then
	#Retrieve selected input database
	wget $dbAddress
	#Extract the database
	dbFile=$(basename $dbAddress)
	gunzip -v $dbFile
	dbFile=$(echo $dbFile | sed 's/\.tgz/\.tar/g')
	tar -xvf $dbFile
elif [[ "$1" == "swissprot" ]]; then
	#Import swissprot database
	"$softsPath"/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
fi
