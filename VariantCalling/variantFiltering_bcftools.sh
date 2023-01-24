#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 8
#$ -N variantFiltering_jobOutput

# script to perform variant filtering after calling
# usage: qsub variantFiltering_bcftools.sh sortedFolderName filterType
# default usage Ex: qsub variantFiltering_bcftools.sh sortedCoordinate_samtoolsHisat2_run1 filteredMapQ
# usage Ex: qsub variantFiltering_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 filteredMapQ
# usage Ex: qsub variantFiltering_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 filteredZS

#Required modules for ND CRC servers
#module load bio

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/"$1
inputsDir=$inputsDir"/variantCallingMerged_"$type

#Retrieve input bam file type
type="$2"

#Make output folder
outFolder=$inputsDir"/variantsFiltered"
mkdir $outFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outFolder"/variantFiltering_summary.txt"
#Name output file of filtering info
outputsFile=$outFolder"/variantFiltering_stats.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

#Check total variants
echo "Total variants from reads with MQ > 60: " > $outputsFile
bcftools filter -i '%QUAL<1001' $inputsDir"/"$type"_calls.vcf.gz" | grep "^scaffold" | wc -l >> $outputsFile

#Include sites with quality > 20 
bcftools filter --threads 8 -i '%QUAL>20' $inputsDir"/"$type"_calls.vcf.gz" -Ob -o $outFolder"/"$type"_calls.flt-qual.bcf"
echo "bcftools filter --threads 8 -i '%QUAL>20' "$inputsDir"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.flt-qual.bcf" >> $inputOutFile
echo "Sites with quality > 20: " >> $outputsFile
bcftools filter --threads 8 -i '%QUAL>20' $inputsDir"/"$type"_calls.vcf.gz" | grep "^scaffold" | wc -l >> $outputsFile

#Include sites with average read depth > 10
bcftools filter --threads 8 -i 'INFO/DP>10' $outFolder"/"$type"_calls.flt-qual.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP.bcf"
echo "bcftools filter --threads 8 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.flt-qual.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf" >> $inputOutFile
echo "Sites with average read depth > 10: " >> $outputsFile
bcftools filter --threads 8 -i 'INFO/DP>10' $outFolder"/"$type"_calls.flt-qual.bcf" | grep "^scaffold" | wc -l >> $outputsFile

#Exclude hetoerozygous sites
bcftools filter --threads 8 -e 'GT="het"' $outFolder"/"$type"_calls.flt-qualDP.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP-homo.bcf"
echo "bcftools filter --threads 8 -e 'GT=\"het\"' "$outFolder"/"$type"_calls.flt-qualDP.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf" >> $inputOutFile
echo "Excluding hetoerozygous sites: " >> $outputsFile
bcftools filter --threads 8 -e 'GT="het"' $outFolder"/"$type"_calls.flt-qualDP.bcf" | grep "^scaffold" | wc -l >> $outputsFile

#Exclude sites homozygous for the reference
bcftools filter --threads 8 -e 'GT="RR"' $outFolder"/"$type"_calls.flt-qualDP-homo.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf"
echo "bcftools filter --threads 8 -e 'GT=\"RR\"' "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf" >> $inputOutFile
echo "Excluding sites homozygous for the reference: " >> $outputsFile
bcftools filter --threads 8 -e 'GT="RR"' $outFolder"/"$type"_calls.flt-qualDP-homo.bcf" | grep "^scaffold" | wc -l >> $outputsFile

#Filter adjacent indels within 2bp
bcftools filter --threads 8 -G 2 $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf" -Ob -o $outFolder"/"$type"_calls.flt-indels.bcf"
echo "bcftools filter --threads 8 -G 2 "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf -Ob -o "$outFolder"/"$type"_calls.flt-indels.bcf" >> $inputOutFile
echo "Excluding adjacent indels within 2bp: " >> $outputsFile
bcftools filter --threads 8 -G 2 $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf" | grep "^scaffold" | wc -l >> $outputsFile

#Filter adjacent SNPs within 2bp
bcftools filter --threads 8 -g 2 $outFolder"/"$type"_calls.flt-indels.bcf" -Ob -o $outFolder"/"$type"_calls.flt-SNPs.bcf"
echo "bcftools filter --threads 8 -g 2 "$outFolder"/"$type"_calls.flt-indels.bcf -Ob -o "$outFolder"/"$type"_calls.flt-SNPs.bcf" >> $inputOutFile
echo "Excluding adjacent SNPs within 2bp: " >> $outputsFile
bcftools filter --threads 8 -g 2 $outFolder"/"$type"_calls.flt-indels.bcf" | grep "^scaffold" | wc -l >> $outputsFile

#Turn on left alignment, normalize indels, and collapse multi allelic sites
bcftools norm --threads 8 -m +any -f $genomeFile" "$outFolder"/"$type"_calls.flt-SNPs.bcf" -Ob -o $outFolder"/"$type"_calls.normCollapse.bcf"
echo "bcftools norm --threads 8 -m +any -f "$genomeFile" "$outFolder"/"$type"_calls.flt-SNPs.bcf -Ob -o "$outFolder"/"$type"_calls.normCollapse.bcf" >> $inputOutFile

#Check final number variants
echo "Turn on left alignment, normalize indels, and collapse multi allelic sites: " >> $outputsFile
bcftools filter -i '%QUAL<1001' $outFolder"/"$type"_calls.normCollapse.bcf" | grep "^scaffold" | wc -l >> $outputsFile
