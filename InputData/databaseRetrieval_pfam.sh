#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_pfam_jobOutput
#Script to retrieve Pfam protein sequence database fasta files
#Note that Pfam databases may be downloaded with wget 
# from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases

#Retreive pfam database storage path
pfamPath=$(grep "pfamDB:" ../InputData/outputPaths.txt | tr -d " " | sed "s/pfamDB://g")
#Move to database directory
cd "$pfamPath"
#Retrieve selected pfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.full.gz
gunzip -v Pfam-A.full.gz