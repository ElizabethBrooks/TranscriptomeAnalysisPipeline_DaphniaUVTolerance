#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_jobOutput
#$ -pe smp 8

#Prepare for mapping
cd ..
dirFlag=0
runNum=0
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	mkdir aligned_hisat2_run"$runNum"
	#Check if the folder already exists
	if [ $? -ne 0 ] ; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for $runNum run of hisat2 alignment..."
	fi
done
mkdir aligned_hisat2_run"$runNum"/out
module load bio
module load bio/hisat2/2.1.0
#Build reference genome
hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_hisat2_run"$runNum"/Daphnia_pulex.allmasked
#Loop through all forward and reverse paired reads and run hisat2 on each pair
# using 8 threads
for f1 in trimmed/*pForward.fq.gz; do
	echo "Sample ${f1:8:${#f1}-23} is being aligned..."
	hisat2 -q -p 8 -x aligned_hisat2/Daphnia_pulex.allmasked -1 $f1 -2 "${f1:0:${#f1}-14}"pReverse.fq.gz -S aligned_hisat2_run"$runNum"/out/"${f1:8:${#f1}-23}".sam --summary-file aligned_hisat2_run"$runNum"/alignedSummary.txt
	#Convert output sam files to bam format for downstream analysis
	echo "Sample ${f1:8:${#f1}-23} is being converted..."
	samtools view -bS aligned_hisat2_run"$runNum"/out/"${f1:8:${#f1}-23}".sam > aligned_hisat2_run"$runNum"/out/"${f1:8:${#f1}-23}".bam
	echo "Sample ${f1:8:${#f1}-23} has been aligned and converted!"
done
