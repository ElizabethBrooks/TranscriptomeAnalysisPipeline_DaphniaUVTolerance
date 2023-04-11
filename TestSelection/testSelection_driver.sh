#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N testSelection_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub testSelection_driver.sh
# usage ex: qsub testSelection_driver.sh

# load necessary modules
module load bio

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# make outputs directory name
outFolder=$inputsPath"/selectionTests"
mkdir $outFolder
#Check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outFolder directory already exsists... please remove before proceeding."
#	exit 1
#fi

# retrieve protein sequences
#bash retrieveFeatures_gffread.sh

# generate Ka and Ks values for protein sequences
bash generateKaKs_musclePal2nalCodeml.sh

