#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantFiltering_jobOutput

# script to perform variant filtering after calling
# usage: qsub variantFiltering_bcftools.sh
# usage Ex: qsub variantFiltering_bcftools.sh

#Required modules for ND CRC servers
module load bio

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsDir=$inputsDir"/variantsMerged_"$type

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

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
echo "Total variants from filtered reads with MQ > 60: " > $outputsFile
bcftools filter -i '%QUAL<1001' $inputsDir"/"$type"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with quality > 20 
bcftools filter --threads 4 -i '%QUAL>20' $inputsDir"/"$type"_calls.norm.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qual.bcf"
echo "bcftools filter --threads 4 -i '%QUAL>20' "$inputsDir"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qual.bcf" >> $inputOutFile
echo "& including sites with quality > 20: " >> $outputsFile
bcftools filter --threads 4 -i '%QUAL>20' $inputsDir"/"$type"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with average read depth > 10
bcftools filter --threads 4 -i 'INFO/DP>10' $outFolder"/"$type"_calls.flt-qual.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP.bcf"
echo "bcftools filter --threads 4 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.flt-qual.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf" >> $inputOutFile
echo "& including sites with average read depth > 10: " >> $outputsFile
bcftools filter --threads 4 -i 'INFO/DP>10' $outFolder"/"$type"_calls.flt-qual.bcf" | grep -v "#" | wc -l >> $outputsFile

#Exclude heterozygous sites
bcftools filter --threads 4 -e 'GT="het"' $outFolder"/"$type"_calls.flt-qualDP.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP-homo.bcf"
echo "bcftools filter --threads 4 -e 'GT=\"het\"' "$outFolder"/"$type"_calls.flt-qualDP.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf" >> $inputOutFile
echo "& excluding heterozygous sites: " >> $outputsFile
bcftools filter --threads 4 -e 'GT="het"' $outFolder"/"$type"_calls.flt-qualDP.bcf" | grep -v "#" | wc -l >> $outputsFile

#Exclude sites homozygous to the reference
bcftools filter --threads 4 -e 'GT="RR"' $outFolder"/"$type"_calls.flt-qualDP-homo.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf"
echo "bcftools filter --threads 4 -e 'GT=\"RR\"' "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf" >> $inputOutFile
echo "& excluding sites homozygous to the reference: " >> $outputsFile
bcftools filter --threads 4 -e 'GT="RR"' $outFolder"/"$type"_calls.flt-qualDP-homo.bcf" | grep -v "#" | wc -l >> $outputsFile

#Filter adjacent indels within 2bp
bcftools filter --threads 4 -G 2 $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf" -Ob -o $outFolder"/"$type"_calls.flt-indels.bcf"
echo "bcftools filter --threads 4 -G 2 "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf -Ob -o "$outFolder"/"$type"_calls.flt-indels.bcf" >> $inputOutFile
echo "& excluding adjacent indels within 2bp: " >> $outputsFile
bcftools filter --threads 4 -G 2 $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf" | grep -v "#" | wc -l >> $outputsFile

#Filter adjacent SNPs within 2bp
bcftools filter --threads 4 -g 2 $outFolder"/"$type"_calls.flt-indels.bcf" -Ob -o $outFolder"/"$type"_calls.flt-SNPs.bcf"
echo "bcftools filter --threads 4 -g 2 "$outFolder"/"$type"_calls.flt-indels.bcf -Ob -o "$outFolder"/"$type"_calls.flt-SNPs.bcf" >> $inputOutFile
echo "& excluding adjacent SNPs within 2bp: " >> $outputsFile
bcftools filter --threads 4 -g 2 $outFolder"/"$type"_calls.flt-indels.bcf" | grep -v "#" | wc -l >> $outputsFile

#Turn on left alignment, normalize indels, and collapse multi allelic sites
bcftools norm --threads 4 -m +any -f $genomeFile $outFolder"/"$type"_calls.flt-SNPs.bcf" -Ob -o $outFolder"/"$type"_calls.normCollapse.bcf"
echo "bcftools norm --threads 4 -m +any -f "$genomeFile" "$outFolder"/"$type"_calls.flt-SNPs.bcf -Ob -o "$outFolder"/"$type"_calls.normCollapse.bcf" >> $inputOutFile
echo "& with left alignment, normalized indels, and collapsed multi allelic sites: " >> $outputsFile
bcftools norm --threads 4 -m +any -f $genomeFile $outFolder"/"$type"_calls.flt-SNPs.bcf" | grep -v "#" | wc -l >> $outputsFile

#Index bcf file
bcftools index --threads 4 $outFolder"/"$type"_calls.normCollapse.bcf"
echo "bcftools index --threads 4 "$outFolder"/"$type"_calls.normCollapse.bcf" >> $inputOutFile

# clean up
rm $outFolder"/"$type"_calls.flt-qual.bcf"
rm $outFolder"/"$type"_calls.flt-qualDP.bcf"
rm $outFolder"/"$type"_calls.flt-qualDP-homo.bcf"
rm $outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf"
rm $outFolder"/"$type"_calls.flt-indels.bcf"
rm $outFolder"/"$type"_calls.flt-SNPs.bcf"
