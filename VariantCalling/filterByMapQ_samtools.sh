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

# retrieve input folder of trimmed data
inputFolder="sortedCoordinate_samtoolsHisat2"

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/"$inputFolder

#Name output file of inputs
inputOutFile=$inputsDir"/mapqFiltering_summary.txt"

#Add version to output file of inputs
samtools --version > $inputOutFile

#Keep only unique read alignments using a mapq score of 60 
#samtools view -@ 4 -bq 60 $inputsDir"/OLYM/accepted_hits.bam" > $inputsDir"/OLYM/filteredMapQ.bam"
#echo "samtools view -@ 8 -bq 60 "$inputsDir"/OLYM/accepted_hits.bam > "$inputsDir"/OLYM/filteredMapQ.bam" >> $inputOutFile
# index bamfile
#samtools index -@ 4 $inputsDir"/OLYM/filteredMapQ.bam"
# clean up
#rm $inputsDir"/OLYM/accepted_hits.bam"

#Keep only unique read alignments using a mapq score of 60 
for f in $inputsDir"/"*"/accepted_hits.bam"; do 
	echo "Processing file $f"
	path=$(dirname $f)
	samtools view -@ 4 -bq 60 $f > $path"/filteredMapQ.bam"
	echo "samtools view -@ 8 -bq 60 "$f" > "$path"/filteredMapQ.bam" >> $inputOutFile
	# index bamfile
	samtools index -@ 4 $path"/filteredMapQ.bam"
done
