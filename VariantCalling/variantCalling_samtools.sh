#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput
#Script to perform variant calling
#Usage: qsub variantCalling_samtools.sh sortedNameFolder analysisTarget
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run3 genome filteredZS
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run1 trimmed_run1E05_assemblyTrinity filteredMapQ
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run1 sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity filteredMapQ

#Required modules for ND CRC servers
module load bio
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
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
inputOutFile="$outputFolder"/counted_summary.txt

#Retrieve input bam file type
type="$3"
#Loop over MapQ filtered bam files
for f in "$inputsDir"/*/"$type".bam; do
	echo "Processing file $f"
	path=$(dirname $f)
	#Calculate the read coverage of positions in the genome
	bcftools mpileup -Ob -o "$path"/"$type"_raw.bcf -f "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta" "$f" 
	#Detect the single nucleotide polymorphisms 
	bcftools call -mv -Ob -o "$path"/"$type"_calls.vcf.gz "$path"/"$type"_raw.bcf 
	#Index vcf file
	bcftools index "$path"/"$type"_calls.vcf.gz
	#Normalize indels
	bcftools norm -f "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta" "$path"/"$type"_calls.vcf.gz -Ob -o "$path"/"$type"_calls.norm.bcf
	#Filter adjacent indels within 5bp
	bcftools filter --IndelGap 5 "$path"/"$type"_calls.norm.bcf -Ob -o "$path"/"$type"_calls.norm.flt-indels.bcf
	#Include sites where FILTER is true
	bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$path"/"$type"_calls.norm.flt-indels.bcf > "$path"/"$type"_filtered_excluded.bcf
done