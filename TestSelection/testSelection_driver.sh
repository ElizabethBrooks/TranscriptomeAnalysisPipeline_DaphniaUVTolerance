#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N testSelection_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub testSelection_driver.sh
# usage ex: qsub testSelectionDriver.sh

# load necessary modules
module load bio

# make outputs directory name
outFolder=$inputsPath"/selectionTests"
mkdir $outFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# retrieve protein sequences
bash retrieveFeatures_gffread.sh

# generate Ka and Ks values for protein sequences
bash generateKaKs_musclePal2nalCodeml.sh

