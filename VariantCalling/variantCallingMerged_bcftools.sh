#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantCallingMerged_jobOutput
#Script to perform variant calling
#Usage: qsub variantCallingMerged_bcftools.sh sortedNameFolder analysisTarget filterType
#Usage Ex: qsub variantCallingMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: qsub variantCallingMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Determine what analysis method was used for the input folder of data
if [[ "$2" == *assemblyTrinity* || "$2" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsDir="$inputsPath"/"$2"/"$1"
	outputsPath="$inputsPath"/"$2"
elif [[ "$2" == *assembly*Trinity* || "$2" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsDir="$inputsPath"/"$2"/"$1"
	outputsPath="$inputsPath"/"$2"
elif [[ "$2" == genome ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsDir="$inputsPath"/"$1"
	outputsPath="$inputsPath"
else
	echo "ERROR: Invalid sorted folder of bam files entered... exiting"
	exit 1
fi

#Set input bam list
inputBamList=../InputData/bamList_Olympics_bcftools.txt
#Set input sample names
#inputSampleList=../InputData/sampleList_Olympics_bcftools.txt

#Make output folder
outFolder="$inputsDir"/variantCallingBcftools_"$3"
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outFolder"/variantCalling_summary.txt

#Select lines associated with input genotype
#genotype=_"$4"_
#grep "$genotype" "$inputBamList" > tmpList_genotype.txt
#Add file type to end of each sample path
type=/"$3".bam
typeTag=$(echo $type | sed "s/\//SLASH/g")
sed -e "s/$/$typeTag/" "$inputBamList" > tmpList.txt
#Add directory to beginning of each sample path
inDir="$inputsDir"/
inDirTag=$(echo $inDir | sed "s/\//SLASH/g")
sed -i -e "s/^/$inDirTag/" tmpList.txt
#Add in slashes
sed -i "s/SLASH/\//g" tmpList.txt

#Retrieve input bam file type
type="$3"

#Output status mesasge
echo "Generating variants for the following input set of bam files: " > "$inputOutFile"
cat tmpList.txt >> "$inputOutFile"

#Calculate the read coverage of positions in the genome
bcftools mpileup --threads 4 -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" -b tmpList.txt
echo bcftools mpileup --threads 4 -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" -b tmpList.txt >> "$inputOutFile"

#Detect the single nucleotide polymorphisms 
bcftools call --threads 4 -mv -Oz -o "$outFolder"/"$type"_calls.vcf.gz "$outFolder"/"$type"_raw.bcf 
echo bcftools call --threads 4 -mv -Oz -o "$outFolder"/"$type"_calls.vcf.gz "$outFolder"/"$type"_raw.bcf >> "$inputOutFile"

#Index vcf file
bcftools index --threads 4 "$outFolder"/"$type"_calls.vcf.gz
echo bcftools index --threads 4 "$outFolder"/"$type"_calls.vcf.gz >> "$inputOutFile"

#Clean up
rm tmpList*.txt