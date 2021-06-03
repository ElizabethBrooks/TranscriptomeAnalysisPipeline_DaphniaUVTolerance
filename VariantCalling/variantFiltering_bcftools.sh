#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantFiltering_jobOutput
#Script to perform variant filtering
#Usage: qsub variantFiltering_bcftools.sh sortedNameFolder analysisTarget variantsFolder
#Usage Ex: qsub variantFiltering_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome variantCallingBcftools_filteredMapQ
#Usage Ex: qsub variantFiltering_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome variantCallingGATK_filteredZS

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
outFolder="$inputsDir"/"$3"
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outFolder"/variantFiltering_summary.txt
#Name output file of filtering info
outputsFile="$outFolder"/variantFiltering_stats.txt

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

#Check total variants
echo "Total variants: " > "$outputsFile"
bcftools filter -i '%QUAL<1001' "$outFolder"/"$type"_calls.vcf.gz | grep "^scaffold" | wc -l >> "$outputsFile"

#Include sites with quality > 20 
bcftools filter --threads 4 -i '%QUAL>20' "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf
echo bcftools filter --threads 4 -i '%QUAL>20' "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf >> "$inputOutFile"
echo "Sites with quality > 20: " >> "$outputsFile"
bcftools filter --threads 4 -i '%QUAL>20' "$outFolder"/"$type"_calls.vcf.gz | grep "^scaffold" | wc -l >> "$outputsFile"

#Include sites with average read depth > 10
bcftools filter --threads 4 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf
echo bcftools filter --threads 4 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.vcf.gz -Ob -o "$outFolder"/"$type"_calls.flt-qualDP.bcf >> "$inputOutFile"
echo "Sites with average read depth > 10: " >> "$outputsFile"
bcftools filter --threads 4 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.vcf.gz | grep "^scaffold" | wc -l >> "$outputsFile"

#Exclude hetoerozygous sites
bcftools filter --threads 4 -e 'GT="het"' "$outFolder"/"$type"_calls.flt-qualDP.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf
echo bcftools filter --threads 4 "-e 'GT=\"het\"'" "$outFolder"/"$type"_calls.flt-qualDP.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf >> "$inputOutFile"
echo "Excluding hetoerozygous sites: " >> "$outputsFile"
bcftools filter --threads 4 -e 'GT="het"' "$outFolder"/"$type"_calls.flt-qualDP.bcf | grep "^scaffold" | wc -l >> "$outputsFile"

#Exclude sites homozygous for the reference
bcftools filter --threads 4 -e 'GT="RR"' "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf
echo bcftools filter --threads 4 "-e 'GT=\"RR\"'" "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf >> "$inputOutFile"
echo "Excluding sites homozygous for the reference: " >> "$outputsFile"
bcftools filter --threads 4 -e 'GT="RR"' "$outFolder"/"$type"_calls.flt-qualDP-homo.bcf | grep "^scaffold" | wc -l >> "$outputsFile"

#Filter adjacent indels within 2bp
bcftools filter --threads 4 -G 2 "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf -Ob -o "$outFolder"/"$type"_calls.flt-indels.bcf
echo bcftools filter --threads 4 -G 2 "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf -Ob -o "$outFolder"/"$type"_calls.flt-indels.bcf >> "$inputOutFile"
echo "Excluding adjacent indels within 2bp: " >> "$outputsFile"
bcftools filter --threads 4 -G 2 "$outFolder"/"$type"_calls.flt-qualDP-homo-dif.bcf | grep "^scaffold" | wc -l >> "$outputsFile"

#Filter adjacent SNPs within 2bp
bcftools filter --threads 4 -g 2 "$outFolder"/"$type"_calls.flt-indels.bcf -Ob -o "$outFolder"/"$type"_calls.flt-SNPs.bcf
echo bcftools filter --threads 4 -g 2 "$outFolder"/"$type"_calls.flt-indels.bcf -Ob -o "$outFolder"/"$type"_calls.flt-SNPs.bcf >> "$inputOutFile"
echo "Excluding adjacent SNPs within 2bp: " >> "$outputsFile"
bcftools filter --threads 4 -g 2 "$outFolder"/"$type"_calls.flt-indels.bcf | grep "^scaffold" | wc -l >> "$outputsFile"

#Turn on left alignment, normalize indels, and collapse multi allelic sites
bcftools norm --threads 4 -m +any -f "$genomeFile" "$outFolder"/"$type"_calls.flt-SNPs.bcf -Ob -o "$outFolder"/"$type"_calls.normCollapse.bcf
echo bcftools norm --threads 4 -m +any -f "$genomeFile" "$outFolder"/"$type"_calls.flt-SNPs.bcf -Ob -o "$outFolder"/"$type"_calls.normCollapse.bcf >> "$inputOutFile"
echo "Turn on left alignment, normalize indels, and collapse multi allelic sites: " >> "$outputsFile"
bcftools filter -i '%QUAL<1001' "$outFolder"/"$type"_calls.normCollapse.bcf | grep "^scaffold" | wc -l >> "$outputsFile"

#Clean up
rm tmpList*.txt