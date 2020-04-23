#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_proteome_jobOutput
#Script to Download the UniProt reference proteomes for all
# organisms below the input taxonomy node ID in compressed FASTA format
#Usage: qsub buildDatabase_proteome.sh taxonID
#Usage ex: qsub buildDatabase_proteome.sh 7215
#Note that the taxonomy ID for the top node of an ogranism
# may be retrieved from https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/

#Load necessary modules
module load bio/blast+
#Set output reference proteome fasta path
outputPath=$(grep "proteomeDB:" inputPaths.txt | tr -d " " | sed "s/proteomeDB://g")
#Make output database directory
outputPath=$(dirname "$outputPath")
mkdir $outputPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	#Error message
	echo "Database directory for taxon $1 already exsists... exiting!"
	exit 1
fi
#Retrieve selected input database
perl ../util/referenceProteomes_byTaxon.pl $1
#Move output fasta to DB folder
#dbFile=$(echo ./*.fasta.gz)
#mv $dbFile $outputPath
#Extract the database
#dbFile=$(basename $dbAddress)
#gunzip -v $dbFile
#Index the database for blastp
#dbFileNoEx=$(echo $dbFile | sed 's/\.gz//')
#makeblastdb -in $dbFileNoEx -dbtype prot
