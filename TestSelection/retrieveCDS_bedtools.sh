#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveCDS_jobOutput

# script to retrieve cds from a referene fasta using a gff
# usage: qsub retrieveCDS_bedtools.sh

#Required modules for ND CRC servers
module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/cds"
mkdir $outFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# set inputs folder
inputsPath=$inputsPath"/variantsMerged_"$type

#Retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# retrieve file name of reference
refTag=$(basename $genomeFile)

#Name output file of inputs
inputOutFile=$outFolder"/retrieveCDS_summary.txt"

#Add version to output file of inputs
bedtools --version > $inputOutFile

# retrieve cds for the reference
bedtools getfasta -name -fi $genomeFile -bed $genomeFeatures -fo $outFolder"/"$refTag"_cds.fa"

# retrieve cds for the consensus
bedtools getfasta -name -fi $inputsPath"/"$type"_consensus.fa" -bed $genomeFeatures -fo $outFolder"/"$type"_consensus_cds.fa"
