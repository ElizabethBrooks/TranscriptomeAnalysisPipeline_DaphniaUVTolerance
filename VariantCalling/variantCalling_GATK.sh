#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 8
#$ -N variantCallingGATK_jobOutput
#Script to perform variant calling
#Usage: qsub variantCallingMerged_bcftools.sh sortedNameFolder analysisTarget filterType
#Usage Ex: qsub variantCallingGATK_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: qsub variantCallingGATK_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS

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
outFolder="$inputsDir"/variantCallingGATK_"$3"
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

while read -r line; do
	#Output status message
	echo "Processing sample: "
	echo "$line"
	
	#Mark duplicates and sort
	java -jar picard.jar MarkDuplicates I="$line" O="$outFolder"/"$type"_mDups.bam M="$outFolder"/"$type"_marked_dup_metrics.txt
	echo java -jar picard.jar MarkDuplicates I="$line" O="$outFolder"/"$type"_mDups.bam M="$outFolder"/"$type"_marked_dup_metrics.txt >> "$inputOutFile"

	#Split reads with N in cigar
	gatk SplitNCigarReads -R "$genomeFile" -I "$outFolder"/"$type"_mDups.bam -O "$outFolder"/"$type"_split.bam
	echo gatk SplitNCigarReads -R "$genomeFile" -I "$outFolder"/"$type"_mDups.bam -O "$outFolder"/"$type"_split.bam >> "$inputOutFile"

	#Generate recalibration table for Base Quality Score Recalibration (BQSR)
	#gatk BaseRecalibrator \
	#	-I "$outFolder"/"$type"_split.bam \
	#	-R "$genomeFile" \
	#	--known-sites sites_of_variation.vcf \
	#	-O "$outFolder"/"$type"_recal_data.table

	#Apply base quality score recalibration
	#gatk ApplyBQSR \
	#	-R "$genomeFile" \
	#	-I "$outFolder"/"$type"_split.bam \
	#	--bqsr-recal-file "$outFolder"/"$type"_recal_data.table \
	#	-O "$outFolder"/"$type"_recal.bam

	#Evaluate and compare base quality score recalibration (BQSR) tables
	#gatk AnalyzeCovariates \
	#    -bqsr "$outFolder"/"$type"_recal_data.table \
	#    -plots "$outFolder"/"$type"_AnalyzeCovariates.pdf

	#Call germline SNPs and indels via local re-assembly of haplotypes
	gatk --java-options "-Xmx4g" HaplotypeCaller  -R "$genomeFile" -I "$outFolder"/"$type"_split.bam -O "$outFolder"/"$type"_hap.g.vcf.gz -ERC GVCF
	echo gatk --java-options "-Xmx4g" HaplotypeCaller  -R "$genomeFile" -I "$outFolder"/"$type"_split.bam -O "$outFolder"/"$type"_hap.g.vcf.gz -ERC GVCF >> "$inputOutFile"
done < tmpList.txt

#Clean up
rm tmpList.txt