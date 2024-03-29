#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N annotation_interproscan_jobOutput
#$ -pe smp 16

# script to perform annotation of protein sequences using interproscan

# usage: qsub annotation_interproscan.sh
# usage Ex: qsub annotation_interproscan.sh

# required modules for ND CRC servers
module load bio

# retrieve protein sequences absolute path
inputsPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")

# retrieve annotation outputs absolute path
outputsPath=$(grep "annotation:" ../InputData/outputPaths.txt | tr -d " " | sed "s/annotation://g")

# set outputs directory name
outputsPath=$outputsPath"/annotated_interproscan"

# make outputs directory
mkdir $outputsPath

# retrieve software path
softwarePath=$(grep "interproscan:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/interproscan://g")

# move to software directory
cd $softwarePath

# run interpro scan on protein sequences
# including GO terms and pathways
./interproscan.sh -cpu 16 -goterms -pa -verbose -vtsv -i $inputsPath -b $outputsPath
