#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantFiltering_jobOutput

# script to perform variant filtering after calling
# usage: qsub variantFiltering_bcftools.sh

# required modules for ND CRC servers
module load bio

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsDir=$inputsDir"/variantsMerged_"$type

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# set output folder
outFolder=$inputsDir

# name output file of inputs
inputOutFile=$outFolder"/variantFiltering_summary.txt"
# name output file of filtering info
outputsFile=$outFolder"/variantFiltering_stats.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

#Check total variants
echo "Total variants from filtered reads with MQ > 60: " > $outputsFile
bcftools view --threads 4 $inputsDir"/"$type"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with quality > 20 
bcftools filter --threads 4 -i '%QUAL>20' $inputsDir"/"$type"_calls.norm.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qual.bcf"
echo "bcftools filter --threads 4 -i '%QUAL>20' "$inputsDir"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qual.bcf" >> $inputOutFile
echo "& including sites with quality > 20: " >> $outputsFile
bcftools view --threads 4 $inputsDir"/"$type"_calls.flt-qual.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with average read depth > 10
bcftools filter --threads 4 -i 'INFO/DP>10' $outFolder"/"$type"_calls.flt-qual.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP.bcf"
echo "bcftools filter --threads 4 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.flt-qual.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf" >> $inputOutFile
echo "& including sites with average read depth > 10: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-qualDP.bcf" | grep -v "#" | wc -l >> $outputsFile

#Exclude heterozygous sites
bcftools filter --threads 4 -e 'GT="het"' $outFolder"/"$type"_calls.flt-qualDP.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP-homo.bcf"
echo "bcftools filter --threads 4 -e 'GT=\"het\"' "$outFolder"/"$type"_calls.flt-qualDP.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf" >> $inputOutFile
echo "& excluding heterozygous sites: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-qualDP-homo.bcf" | grep -v "#" | wc -l >> $outputsFile

# remove variants with uncalled genotypes
bcftools filter --threads 4 -i 'INFO/AN=8' $outFolder"/"$type"_calls.flt-qualDP-homo.bcf" -Ob -o $outFolder"/"$type"_calls.flt-AN.bcf"
echo "bcftools filter --threads 4 -i 'INFO/AN=8' "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf -Ob -o "$outFolder"/"$type"_calls.flt-AN.bcf" >> $inputOutFile
echo "& excluding uncalled genotypes: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-AN.bcf" | grep -v "#" | wc -l >> $outputsFile

# keep variants that have a MQSB value
# Mann-Whitney U test of Mapping Quality vs Strand Bias
# see MQB in https://www.reneshbedre.com/blog/vcf-fields.html
bcftools filter --threads 4 -i 'INFO/MQSB>0.05' $outFolder"/"$type"_calls.flt-AN.bcf" -Ob -o $outFolder"/"$type"_calls.flt-MQSB.bcf"
echo "bcftools filter --threads 4 -i 'INFO/MQSB>0.05' "$outFolder"/"$type"_calls.flt-AN.bcf -Ob -o "$outFolder"/"$type"_calls.flt-MQSB.bcf" >> $inputOutFile
echo "& including variants with MQSB > 0.05: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-MQSB.bcf" | grep -v "#" | wc -l >> $outputsFile

#Turn on left alignment, normalize indels, and collapse multi allelic sites
bcftools norm --threads 4 -m +any -f $genomeFile $outFolder"/"$type"_calls.flt-MQSB.bcf" -Ob -o $outFolder"/"$type"_calls.flt-norm.bcf"
echo "bcftools norm --threads 4 -m +any -f "$genomeFile" "$outFolder"/"$type"_calls.flt-MQSB.bcf -Ob -o "$outFolder"/"$type"_calls.flt-norm.bcf" >> $inputOutFile
echo "& with left alignment, normalized indels, and collapsed multi allelic sites: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Index bcf file
bcftools index --threads 4 $outFolder"/"$type"_calls.flt-norm.bcf"
echo "bcftools index --threads 4 "$outFolder"/"$type"_calls.flt-norm.bcf" >> $inputOutFile

# clean up
#rm $outFolder"/"$type"_calls.flt-qual.bcf"
#rm $outFolder"/"$type"_calls.flt-qualDP.bcf"
#rm $outFolder"/"$type"_calls.flt-qualDP-homo.bcf"
#rm $outFolder"/"$type"_calls.flt-AN.bcf"
#rm $outFolder"/"$type"_calls.flt-MQSB.bcf"
