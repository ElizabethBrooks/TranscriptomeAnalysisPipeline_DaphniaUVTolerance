#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N variantPrepGATK_jobOutput
#Script to perform variant prep
#Usage: qsub variantPrep_GATK.sh sortedNameFolder analysisTarget filterType
#Usage Ex: qsub variantPrep_GATK.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: qsub variantPrep_GATK.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS

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
elif [[ "$2" == *assembly*Trinity* || "$2" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsDir="$inputsPath"/"$2"/"$1"
elif [[ "$2" == genome ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsDir="$inputsPath"/"$1"
else
	echo "ERROR: Invalid sorted folder of bam files entered... exiting"
	exit 1
fi
#Retrieve outputs absolute path
outputsPath=$(grep "scratch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/scratch://g")
outputsDir="$outputsPath"/"$1"_"$2"
mkdir "$outputsDir"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsDir directory already exsists... please remove before proceeding."
	exit 1
fi

#Set input bam list
inputBamList=../InputData/fileList_Olympics.txt

#Name output file of inputs
inputOutFile="$outFolder"/variantCalling_summary.txt

#Create first outputs directory
outFolder="$outputsDir"/variantCallingGATK_"$3"
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Make second output folder
outFolder="$outFolder"/variantsPreped
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outFolder"/variantPrep_summary.txt

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
echo "Generating variants for the following input set of bam files: " >> "$inputOutFile"
cat tmpList.txt >> "$inputOutFile"

while read -r line; do
	#Clean up sample tag
	tag=$(dirname $line)
	tag=$(basename $tag)

	#Output status message
	echo "Processing sample: "
	echo "$tag"

	#Mark duplicates and sort
	picard MarkDuplicates I="$line" O="$outFolder"/"$tag"_mDups.bam M="$outFolder"/"$tag"_marked_dup_metrics.txt
	echo picard MarkDuplicates I="$line" O="$outFolder"/"$tag"_mDups.bam M="$outFolder"/"$tag"_marked_dup_metrics.txt >> "$inputOutFile"

	#Split reads with N in cigar
	gatk SplitNCigarReads -R "$genomeFile" -I "$outFolder"/"$tag"_mDups.bam -O "$outFolder"/"$tag"_split.bam
	echo gatk SplitNCigarReads -R "$genomeFile" -I "$outFolder"/"$tag"_mDups.bam -O "$outFolder"/"$tag"_split.bam >> "$inputOutFile"

	#Generate recalibration table for Base Quality Score Recalibration (BQSR)
	#gatk BaseRecalibrator -I "$outFolder"/"$tag"_split.bam -R "$genomeFile" --known-sites sites_of_variation.vcf -O "$outFolder"/"$tag"_recal_data.table

	#Apply base quality score recalibration
	#gatk ApplyBQSR -R "$genomeFile" -I "$outFolder"/"$tag"_split.bam --bqsr-recal-file "$outFolder"/"$tag"_recal_data.table -O "$outFolder"/"$tag"_recal.bam

	#Evaluate and compare base quality score recalibration (BQSR) tables
	#gatk AnalyzeCovariates -bqsr "$outFolder"/"$tag"_recal_data.table -plots "$outFolder"/"$tag"_AnalyzeCovariates.pdf

	#Add read group info
	picard AddOrReplaceReadGroups I="$outFolder"/"$tag"_split.bam O="$outFolder"/"$tag"_RG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$tag"
	echo picard AddOrReplaceReadGroups I="$outFolder"/"$tag"_split.bam O="$outFolder"/"$tag"_RG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$tag" >> "$inputOutFile"
done < tmpList.txt

#Clean up
rm tmpList.txt