#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N filterMapQ_jobOutput

# script to perform read quaity filtering of coordinate sorted bam files
# usage: qsub filterByMapQ_samtools.sh

#Required modules for ND CRC servers
module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_samtoolsBcftools"
outputsPath=$inputsPath

#Name output file of inputs
inputOutFile=$inputsDir"/mapqFiltering_summary.txt"

#Add version to output file of inputs
samtools --version > $inputOutFile

#Keep only unique read alignments using a mapq score of 60 
samtools view -@ 4 -bq 60 $inputsDir"/OLYM/accepted_hits.bam" > $inputsDir"/OLYM/accepted_hits_readGroups.bam"
echo "samtools view -@ 8 -bq 60 "$inputsDir"/OLYM/accepted_hits.bam > "$inputsDir"/OLYM/accepted_hits_readGroups.bam" >> $inputOutFile
# index bamfile
samtools index -@ 4 $inputsDir"/OLYM/accepted_hits_readGroups.bam"
# clean up
rm $inputsDir"/OLYM/accepted_hits.bam"

#Keep only unique read alignments using a mapq score of 60 
#for f in $inputsDir"/"*"/accepted_hits_RG.bam"; do 
#	echo "Processing file $f"
#	path=$(dirname $f)
#	samtools view -@ 4 -bq 60 $f > $path"/accepted_hits_readGroups.bam"
#	echo "samtools view -@ 8 -bq 60 "$f" > "$path"/accepted_hits_readGroups.bam" >> $inputOutFile
	# index bamfile
#	samtools index -@ 4 $path"/accepted_hits_readGroups.bam"
	# clean up
#	rm $f
#done
