#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N database_uniprot_jobOutput
#Script to retrieve Uniprot protein sequence database fasta files
#Note that uniprot databases may be downloaded from the UniprotKB search page, or wget
#ex: wget https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:arthropod
#or ex: wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz

#Retreive uniprot database storage path
uniprotPath=$(grep "uniprotDB:" inputPaths.txt | tr -d " " | sed "s/uniprotDB://g")
uniprotPath=$(dirname "$uniprotPath")
#Move to database directory
cd "$uniprotPath"
#Retrieve selected uniprot database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
gunzip -v uniref100.fasta.gz