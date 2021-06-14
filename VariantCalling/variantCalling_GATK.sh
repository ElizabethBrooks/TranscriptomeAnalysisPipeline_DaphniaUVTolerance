#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N variantCallingGATK_jobOutput
#Script to perform variant calling
#Usage: qsub variantCalling_GATK.sh sortedNameFolder analysisTarget
#Usage Ex: qsub variantCalling_GATK.sh sortedCoordinate_samtoolsHisat2_run3 variantCallingGATK_filteredMapQ

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve inputs absolute path
inputsPath=$(grep "scratch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/scratch://g")
inputsDir="$inputsPath"/"$1"/"$2"/variantsPreped

#Set input file list
inputBamList=../InputData/fileList_Olympics.txt

#Make output folder
outFolder=$(dirname $inputsDir)
outFolder="$outFolder"/variantsCalled
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outFolder"/variantCalling_summary.txt

#Output status mesasge
echo "Generating variants for the following input set of bam files: " > "$inputOutFile"
cat tmpList.txt >> "$inputOutFile"

while read -r line; do
	#Output status message
	echo "Processing sample: "
	echo "$line"_RG.bam

	#Index input bam
	samtools index "$line"_RG.bam

	#Call germline SNPs and indels via local re-assembly of haplotypes
	gatk --java-options "-Xmx4g" HaplotypeCaller  -R "$genomeFile" -I "$inputsDir"/"$line"_RG.bam -O "$outFolder"/"$line"_hap.g.vcf.gz -ERC GVCF
	echo gatk --java-options "-Xmx4g" HaplotypeCaller  -R "$genomeFile" -I "$inputsDir"/"$line"_RG.bam -O "$outFolder"/"$line"_hap.g.vcf.gz -ERC GVCF >> "$inputOutFile"
done < tmpList.txt

#Clean up
rm tmpList.txt