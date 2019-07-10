#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
module load bio/hisat2/2.1.0
#Prepare for mapping
cd ..
dirFlag=0
runNum=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve folders to analyze from the input arguments
for f1 in "$@"; do
	#Make a new directory for each alignment run
	while [ $dirFlag -eq 0 ]; do
		mkdir aligned_hisat2_run"$runNum"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of hisat2 alignment on $f1 data..."
		fi
	done
	#Build reference genome
	hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_hisat2_run"$runNum"/Daphnia_pulex.allmasked
	#Loop through all forward and reverse paired reads and run hisat2 on each pair
	# using 8 threads and samtools to convert output sam files to bam
	for f2 in "$f1"/*pForward.fq.gz; do
		echo "Samples $f2 and ${f2:13:${#f2}-28} are being aligned and converted..."
		hisat2 -p 8 -q -x aligned_hisat2_run"$runNum"/Daphnia_pulex.allmasked -1 $f2 -2 "${f2:13:${#f2}-28}"pReverse.fq.gz -S aligned_hisat2_run"$runNum"/out/"${f2:13:${#f2}-28}".sam --summary-file aligned_hisat2_run"$runNum"/alignedSummary.txt | samtools view -@ 8 -bS aligned_hisat2_run"$runNum"/out/"${f2:13:${#f2}-28}".sam > aligned_hisat2_run"$runNum"/out/"${f2:13:${#f2}-28}".bam
		#hisat2 -p 8 -q -x aligned_hisat2_run"$runNum"/Daphnia_pulex.allmasked -1 $f2 -2 "${f2:13:${#f2}-28}"pReverse.fq.gz -S aligned_hisat2_run"$runNum"/out/"${f2:13:${#f2}-28}".sam --summary-file aligned_hisat2_run"$runNum"/alignedSummary.txt
		#Convert output sam files to bam format for downstream analysis
		#echo "Samples $f2 and ${f2:13:${#f2}-28} are being converted..."
		#samtools view -@ 8 -bS aligned_hisat2_run"$runNum"/out/"${f2:13:${#f2}-28}".sam > aligned_hisat2_run"$runNum"/out/"${f2:13:${#f2}-28}".bam
		echo "Samples $f2 and ${f2:13:${#f2}-28} have been aligned and converted!"
	done
done
