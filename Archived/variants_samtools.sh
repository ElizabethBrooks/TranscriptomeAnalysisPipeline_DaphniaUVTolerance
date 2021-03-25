#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N variants_samtools_jobOutput
#Script to perform variant calling with bcftools
#Usage: qsub variants_samtools.sh sortedFolder
#Usage Ex: qsub variants_samtools.sh sortedCoordinate_samtoolsTophat2_run1

#Required modules for ND CRC servers
module load bio/2.0
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
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_variants_summary.txt
#Loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/"$1"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Perform variant calling using Samtools bcftools
	#Also normalize the vcf file
	echo "Sample $curSampleNoPath variant are being called..."
	bcftools mpileup -Ou -f "$genomeFile" "$curAlignedSample" | bcftools call -Ou -mv | bcftools norm -f "$genomeFile" Oz -o "$outputFolder"/"$curSampleNoPath".vcf.gz
	#Finally, index the vcf file
	tabix "$outputFolder"/"$curSampleNoPath".vcf.gz
	#Add run inputs to output summary file
	echo "$curSampleNoPath" >> "$inputOutFile"
	echo "bcftools mpileup -Ou -f" "$genomeFile" "$curAlignedSample" "| bcftools call -Ou -mv | bcftools norm -f" "$genomeFile" "Oz -o" "$outputFolder""/""$curSampleNoPath"".vcf.gz" >> "$inputOutFile"
	echo "Sample $curSampleNoPath variants have been called!"
done
#Merge vcf files to genreate genotypes for all individuals at all of the unique positions present across the files
#bcftools merge 1_vcf.gz 2_vcf.gz --threads 16 --missing-to-ref --merge both -O z -o Full_merged.vcf_new.gz
#Create a consesnsus sequence by applying vcf variants to a reference fasta file
#bcftools consensus -f reference.fasta -o out.fa Full_merged.vcf_new.gz
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"
