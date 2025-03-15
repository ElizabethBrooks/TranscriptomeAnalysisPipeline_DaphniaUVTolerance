#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 8
#$ -N variantCalling_jobOutput
#Script to perform variant calling
#Usage: qsub variantCalling_bcftools.sh sortedNameFolder analysisTarget
#Usage Ex: qsub variantCalling_bcftools.sh sortedCoordinate_samtoolsHisat2_run1 genome filteredMapQ
#Usage Ex: qsub variantCalling_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: qsub variantCalling_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS

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
#Name output file of inputs
inputOutFile="$inputsDir"/variantCalling_summary.txt

#Add version to output file of inputs
bcftools --version > "$inputOutFile"

#Retrieve input bam file type
type="$3"
#Loop over MapQ filtered bam files
for f in "$inputsDir"/*/"$type".bam; do
	echo "Processing file $f"
	path=$(dirname $f)
	#Calculate the read coverage of positions in the genome
	bcftools mpileup --threads 8 -d 8000 -Q 20 -Ob -o "$path"/"$type"_raw.bcf -f "$genomeFile" "$f" 
	echo bcftools mpileup --threads 8 -d 8000 -Q 20 -Ob -o "$path"/"$type"_raw.bcf -f "$genomeFile" "$f" >> "$inputOutFile"
	#Detect the single nucleotide polymorphisms 
	bcftools call --threads 8 -mv -Oz -o "$path"/"$type"_calls.vcf.gz "$path"/"$type"_raw.bcf 
	echo bcftools call --threads 8 -mv -Oz -o "$path"/"$type"_calls.vcf.gz "$path"/"$type"_raw.bcf >> "$inputOutFile"
	#Index vcf file
	bcftools index --threads 8 "$path"/"$type"_calls.vcf.gz
	echo bcftools index --threads 8 "$path"/"$type"_calls.vcf.gz >> "$inputOutFile"
	#Normalize indels
	bcftools norm --threads 8 -f "$genomeFile" "$path"/"$type"_calls.vcf.gz -Ob -o "$path"/"$type"_calls.norm.bcf
	echo bcftools norm --threads 8 -f "$genomeFile" "$path"/"$type"_calls.vcf.gz -Ob -o "$path"/"$type"_calls.norm.bcf >> "$inputOutFile"
	#Filter adjacent indels within 5bp
	bcftools filter --threads 8 --IndelGap 5 "$path"/"$type"_calls.norm.bcf -Ob -o "$path"/"$type"_calls.norm.flt-indels.bcf
	echo bcftools filter --threads 8 --IndelGap 5 "$path"/"$type"_calls.norm.bcf -Ob -o "$path"/"$type"_calls.norm.flt-indels.bcf >> "$inputOutFile"
	#Include sites where FILTER is true
	#bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$path"/"$type"_calls.norm.flt-indels.bcf > "$path"/"$type"_filtered.bcf
	#echo bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$path"/"$type"_calls.norm.flt-indels.bcf ">" "$path"/"$type"_filtered.bcf >> "$inputOutFile"
done