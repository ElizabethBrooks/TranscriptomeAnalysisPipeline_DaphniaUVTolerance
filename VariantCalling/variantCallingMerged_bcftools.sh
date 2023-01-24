#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 8
#$ -N variantCallingMerged_jobOutput

# script to perform variant calling of mapq filtered bam files before variant filtering
# usage: qsub variantCallingMerged_bcftools.sh sortedFolderName filterType
# default usage Ex: qsub variantCallingMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run1 filteredMapQ
# usage Ex: qsub variantCallingMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 filteredMapQ
# usage Ex: qsub variantCallingMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 filteredZS

#Required modules for ND CRC servers
#module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/"$1
outputsPath="$inputsPath"

#Set input bam list
inputBamList=../InputData/fileList_Olympics.txt
#Set input sample names
#inputSampleList=../InputData/sampleList_Olympics_bcftools.txt

#Retrieve input bam file type
type="$2"

#Make output folder
outFolder=$inputsDir"/variantCallingMerged_"$type
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outFolder"/variantCalling_summary.txt"

#Add version to output file of inputs
bcftools --version > "$inputOutFile"

#Select lines associated with input genotype
#genotype=_"$4"_
#grep "$genotype" "$inputBamList" > tmpList_genotype.txt
#Add file type to end of each sample path
type="/"$3".bam"
typeTag=$(echo $type | sed "s/\//SLASH/g")
sed -e "s/$/$typeTag/" $inputBamList > "tmpList.txt"
#Add directory to beginning of each sample path
inDir="$inputsDir"/
inDirTag=$(echo $inDir | sed "s/\//SLASH/g")
sed -i -e "s/^/$inDirTag/" "tmpList.txt"
#Add in slashes
sed -i "s/SLASH/\//g" "tmpList.txt"

#Output status mesasge
echo "Generating variants for the following input set of bam files: " >> $inputOutFile
cat tmpList.txt >> "$inputOutFile"

#Calculate the read coverage of positions in the genome
#bcftools mpileup --threads 8 -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" -b tmpList.txt
bcftools mpileup --threads 8 -d 8000 -Q 20 -Ob -o $outFolder"/"$type"_raw.bcf" -f $genomeFile -b "tmpList.txt"
echo "bcftools mpileup --threads 8 -d 8000 -Q 20 -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" "$f >> $inputOutFile

#Detect the single nucleotide polymorphisms 
bcftools call --threads 8 -mv -Oz -o $outFolder"/"$type"_calls.vcf.gz" $outFolder"/"$type"_raw.bcf" 
echo "bcftools call --threads 8 -mv -Oz -o "$outFolder"/"$type"_calls.vcf.gz "$outFolder"/"$type"_raw.bcf" >> $inputOutFile

#Index vcf file
bcftools index --threads 8 $outFolder"/"$type"_calls.vcf.gz"
echo "bcftools index --threads 8 "$outFolder"/"$type"_calls.vcf.gz" >> $inputOutFile

#Include sites where FILTER is true
#bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$outFolder"/"$type"_calls.norm.flt-indels.bcf > "$outFolder"/"$type"_filtered.bcf
#echo bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$outFolder"/"$type"_calls.norm.flt-indels.bcf ">" "$outFolder"/"$type"_filtered.bcf >> "$inputOutFile"

#Clean up
#rm "tmpList.txt"