#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantCallingMerged_jobOutput

# script to perform variant calling of mapq filtered bam files before variant filtering
# usage: qsub variantCallingMerged_bcftools.sh
# usage Ex: qsub variantCallingMerged_bcftools.sh

# load required modules for ND CRC servers
#module load bio

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# set input bam file type
type="filteredMapQ"

# set input bam list
inputBamList="../InputData/fileList_Olympics_genotypes.txt"

#Make output folder
outFolder=$inputsDir"/variantsMerged_"$type
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
typeTag="SLASH"$type".bam"
sed -e "s/$/$typeTag/" $inputBamList > $outFolder"/tmpList.txt"
#Add directory to beginning of each sample path
inDirTag=$(echo $inputsDir"/" | sed "s/\//SLASH/g")
sed -i -e "s/^/$inDirTag/" $outFolder"/tmpList.txt"
#Add in slashes
sed -i "s/SLASH/\//g" $outFolder"/tmpList.txt"

#Output status mesasge
echo "Generating variants for the following input set of bam files: " >> $inputOutFile
cat $outFolder"/tmpList.txt" >> "$inputOutFile"

#Calculate the read coverage of positions in the genome
#bcftools mpileup --threads 4 -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" -b $outFolder"/tmpList.txt"
bcftools mpileup --threads 4 -d 8000 -Q 20 -Ob -o $outFolder"/"$type"_raw.bcf" -f $genomeFile -b $outFolder"/tmpList.txt"
echo "bcftools mpileup --threads 4 -d 8000 -Q 20 -Ob -o "$outFolder"/"$type"_raw.bcf -f "$genomeFile" "$f >> $inputOutFile

#Detect the single nucleotide polymorphisms 
bcftools call --threads 4 -mv -Oz -o $outFolder"/"$type"_calls.vcf.gz" $outFolder"/"$type"_raw.bcf" 
echo "bcftools call --threads 4 -mv -Oz -o "$outFolder"/"$type"_calls.vcf.gz "$outFolder"/"$type"_raw.bcf" >> $inputOutFile

#Index vcf file
bcftools index --threads 4 $outFolder"/"$type"_calls.vcf.gz"
echo "bcftools index --threads 4 "$outFolder"/"$type"_calls.vcf.gz" >> $inputOutFile

#Normalize indels
bcftools norm --threads 4 -f $genomeFile $outFolder"/"$type"_calls.vcf.gz" -Ob -o $outFolder"/"$type"_calls.norm.bcf"
echo "bcftools norm --threads 4 -f "$genomeFile" "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.norm.bcf" >> "$inputOutFile"

#Filter adjacent indels within 5bp
#bcftools filter --threads 4 --IndelGap 5 "$outFolder"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.norm.flt-indels.bcf
#echo bcftools filter --threads 4 --IndelGap 5 "$outFolder"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.norm.flt-indels.bcf >> "$inputOutFile"

#Include sites where FILTER is true
#bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$outFolder"/"$type"_calls.norm.flt-indels.bcf > "$outFolder"/"$type"_filtered.bcf
#echo bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$outFolder"/"$type"_calls.norm.flt-indels.bcf ">" "$outFolder"/"$type"_filtered.bcf >> "$inputOutFile"

#Clean up
rm $outFolder"/tmpList.txt"