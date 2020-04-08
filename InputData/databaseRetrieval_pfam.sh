#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_pfam_jobOutput
#Script to retrieve Pfam protein sequence database fasta files
#Usage: qsub databaseRetrieval_pfam.sh 
#Note that Pfam databases may be downloaded with wget 
# from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases

#Retreive pfam database storage path
pfamPath=$(grep "pfamDB:" inputPaths.txt | tr -d " " | sed "s/pfamDB://g")
pfamPath=$(dirname "$pfamPath")
#Move to database directory
cd "$pfamPath"
#Retrieve selected pfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz
gunzip -v Pfam-A.hmm.gz