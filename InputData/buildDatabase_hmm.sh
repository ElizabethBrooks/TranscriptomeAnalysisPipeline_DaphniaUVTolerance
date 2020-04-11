#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_hmm_jobOutput
#Script to retrieve Pfam protein sequence database fasta files
#Usage: qsub databaseRetrieval_hmm.sh 
#Usage ex: qsub databaseRetrieval_hmm.sh 
#Note that Pfam databases may be downloaded with wget 
# from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases

#Load necessary modules for ND CRC servers
module load bio/hmmer
#Retreive pfam database storage path
outputPath=$(grep "pfamDB:" inputPaths.txt | tr -d " " | sed "s/pfamDB://g")
#Make output database directory
outputPath=$(dirname "$outputPath")
mkdir $outputPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	#Error message
	echo "Database directory for $1 files already exsists... exiting!"
	exit 1
fi
#Move to database directory
cd $outputPath
#Retrieve selected pfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz
gunzip -v Pfam-A.hmm.gz
#Use hmmpress to prepare a HMM database for hmmscan
hmmpress Pfam-A.hmm