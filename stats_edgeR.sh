#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2
#Prepare for analysis
cd ..
dirFlag=0
runNum=0
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	mkdir stats_edgeR_run"$runNum"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for $runNum run of edgeR stats analysis..."
	fi
done
mkdir stats_edgeR_run"$runNum"/sorted
genomeFile=TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/PA42.3.0.annotation.18440.gff
#Loop through all forward and reverse paired reads and store the file locations in arrays
for f1 in aligned_tophat2/out/*; do
	echo "Sample ${f1:20:${#f1}-0} is being analyzed..."
	#Run samtools to prepare mapped reads for counting
	# using 8 threads
	#samtools sort -o stats_edgeR_run"$runNum"/sorted/${f1:20:${#f1}-0}/accepted_hits.sorted.bam -T /tmp/${f1:20:${#f1}-0}/accepted_hits.sorted.bam -@ 8 $f1/accepted_hits.bam
	#Run htseq-count to prepare sorted reads for stats analysis in edgeR
	#htseq-count -s no -m union -t gene -i trID $f1/accepted_hits.sorted.bam -i "$genomeFile" > stats_edgeR_run"$runNum"/${f1:20:${#f1}-0}.counts
done
