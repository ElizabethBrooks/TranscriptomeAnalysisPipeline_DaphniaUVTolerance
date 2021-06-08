#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantCallingGATK_jobOutput
#Script to perform variant calling
#Usage: qsub variantCombo_GATK.sh sortedNameFolder analysisTarget filterType
#Usage Ex: qsub variantCombo_GATK.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: qsub variantCombo_GATK.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve genome features absolute path for alignment
scratchPath=$(grep "scratch" ../InputData/outputPaths.txt | tr -d " " | sed "s/scratch://g")
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
inputSampleList=../InputData/sampleList_Olympics_bcftools.txt

#Make output folder
outFolder="$inputsDir"/variantCallingGATK_"$3"
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outFolder"/variantCombo_summary.txt

#Select lines associated with input genotype
#genotype=_"$4"_
#grep "$genotype" "$inputBamList" > tmpList_genotype.txt
#Add file type to end of each sample path
type=/"$3".bam
typeTag=$(echo $type | sed "s/\//SLASH/g")
sed -e "s/$/$typeTag/" "$inputBamList" > tmpList1.txt
#Add directory to beginning of each sample path
inDir="$inputsDir"/
inDirTag=$(echo $inDir | sed "s/\//SLASH/g")
sed -i -e "s/^/$inDirTag/" tmpList1.txt
#Add in slashes
sed -i "s/SLASH/\//g" tmpList1.txt

#Add sample names to each path line
cut -f1 "$inputSampleList" > tmpList2.txt
paste -d tmpList1.txt tmpList2.txt > tmpList_full.txt

#Output status mesasge
echo "Generating variants for the following input set of bam files: " > "$inputOutFile"
cat tmpList_full.txt >> "$inputOutFile"

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path "$outFolder" -L scaffold_1-scaffold_496 --sample-name-map tmpList_full.txt --tmp-dir "$scratchPath" --reader-threads 4
echo gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path "$outFolder" -L scaffold_1-scaffold_496 --sample-name-map tmpList_full.txt --tmp-dir "$scratchPath" --reader-threads 4 >> "$inputOutFile"

#Clean up
rm tmpList*.txt
