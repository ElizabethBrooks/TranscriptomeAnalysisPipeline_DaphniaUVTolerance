#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N variants_GATK_jobOutput
#Script to perform variant calling with gatk:
# SplitNTrim > BQSR > HaplotypeCaller
#Bam files need to have duplicates marked using samtools
#Usage: qsub variants_GATK.sh sortedFolder
#Usage Ex: qsub variants_GATK.sh sortedCoordinate_samtoolsTophat2_run1

#Required modules for ND CRC servers
module load bio/samtools
module load bio/gatk
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != sortedCoordinate* ]]; then
	echo "ERROR: The "$1" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$1" == *"Hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$1" == *"Tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The sorted "$1" folder of bam files were not found... exiting"
	exit 1
fi
#TODO: check input folder requirements
#Retrieve aligned reads input absolute path
inputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve variant calling outputs absolute path
outputsPath=$(grep "variantCalling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/variantCalling://g")
#Create output directory
outputFolder="$outputsPath"/"$1"_variants
mkdir "$outputFolder"
#Move to outputs directory
cd "$outputFolder"
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
#Loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/"$1"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Mark duplicates using smatools fixmate
	samtools fixmate -m "$f1" "$outputFolder"/"$curSampleNoPath"_fixed.bam
	#Remove duplicate reads wtih samtools markdup
	samtools markdup -r "$outputFolder"/"$curSampleNoPath"_fixed.bam "$outputFolder"/"$curSampleNoPath"_marked.bam
	rm "$outputFolder"/"$curSampleNoPath"_fixed.bam
	#Index the marked file
	samtools index "$outputFolder"/"$curSampleNoPath"_marked.bam "$outputFolder"/"$curSampleNoPath"_indexed.bam
	rm "$outputFolder"/"$curSampleNoPath"_marked.bam
	#Perform variant calling using gatk
	echo "Sample $curSampleNoPath variant are being called..."
	#Splits reads into exon segments (getting rid of Ns but maintaining grouping information)
	# and hard­clip any sequences overhanging into the intronic regions
	gatk SplitNCigarReads -R "$genomeFile" -I "$outputFolder"/"$curSampleNoPath"_indexed.bam -O "$outputFolder"/"$curSampleNoPath"_split.bam
	rm "$outputFolder"/"$curSampleNoPath"_indexed.bam
	#Correct for systematic bias that affect the assignment of base quality scores by the sequencer
	gatk ApplyBQSR -R "$genomeFile" -I "$outputFolder"/"$curSampleNoPath"_split.bam --bqsr-recal-file "$outputFolder"/"$curSampleNoPath"_recalibration.table -O "$outputFolder"/"$curSampleNoPath"_recal.bam
	rm "$outputFolder"/"$curSampleNoPath"_split.bam
	#Finally, run HaplotypeCaller in GVCF mode so multiple samples may be added in a scalable way
	gatk HaplotypeCaller -R "$genomeFile" -I "$outputFolder"/"$curSampleNoPath"_recal.bam -O "$outputFolder"/variants.g.vcf -ERC GVCF
	rm "$outputFolder"/"$curSampleNoPath"_recal.bam
done
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"
