#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N variantCallingMerged_jobOutput
#Script to perform variant calling
#Usage: qsub variantCallingMerged_bcftools.sh sortedNameFolder analysisTarget
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

#Make output folder
outFolder="$inputsDir"/variantCalling_"$3"
mkdir "$outFolder"
#Name output file of inputs
inputOutFile="$outFolder"/variantCalling_summary.txt

#Add file type to end of each sample path
typeTag=/"$3".bam
sed -e 's/$/$typeTag/' "$inputBamList" > tmpList.txt
#Add directory to beginning of each sample path
inDirTag="$inputsDir"/
sed -i -e 's/^/$inDirTag/' tmpList.txt

#Output status mesasge
echo "Generating variants for input set of bam files..."

#Calculate the read coverage of positions in the genome
bcftools mpileup -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" -b tmpList.txt
echo bcftools mpileup -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" -b tmpList.txt >> "$inputOutFile"
#Detect the single nucleotide polymorphisms 
bcftools call -mv -Ob -o "$outFolder"/"$type"_calls.vcf.gz "$outFolder"/"$type"_raw.bcf 
echo bcftools call -mv -Ob -o "$outFolder"/"$type"_calls.vcf.gz "$outFolder"/"$type"_raw.bcf >> "$inputOutFile"
#Index vcf file
bcftools index "$outFolder"/"$type"_calls.vcf.gz
echo bcftools index "$outFolder"/"$type"_calls.vcf.gz >> "$inputOutFile"
#Normalize indels
bcftools norm -f "$genomeFile" "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.norm.bcf
echo bcftools norm -f "$genomeFile" "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.norm.bcf >> "$inputOutFile"
#Filter adjacent indels within 5bp
bcftools filter --IndelGap 5 "$outFolder"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.norm.flt-indels.bcf
echo bcftools filter --IndelGap 5 "$outFolder"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.norm.flt-indels.bcf >> "$inputOutFile"
#Include sites where FILTER is true
bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$outFolder"/"$type"_calls.norm.flt-indels.bcf > "$outFolder"/"$type"_filtered.bcf
echo bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$outFolder"/"$type"_calls.norm.flt-indels.bcf ">" "$outFolder"/"$type"_filtered.bcf >> "$inputOutFile"

#Clean up
rm tmpList.txt