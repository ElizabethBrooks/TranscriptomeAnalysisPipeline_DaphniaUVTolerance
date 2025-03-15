#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantComboGATK_jobOutput
#Script to combine variants across samples
#Usage: qsub variantCombo_GATK.sh sortedNameFolder analysisTarget filterType numScaffolds
#Usage Ex: qsub variantCombo_GATK.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ 496
#Usage Ex: qsub variantCombo_GATK.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS 496

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve genome features absolute path for alignment
scratchPath=$(grep "scratch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/scratch://g")
#Retrieve inputs absolute path
inputsPath=$(grep "scratch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/scratch://g")
inputsDir="$inputsPath"/"$1"/"$2"

#Set input file and sample lists
inputsDir="$inputsDir"/variantCallingGATK_"$3"
inputFileList=../InputData/fileList_Olympics.txt
inputSampleList=../InputData/sampleList_Olympics.txt

#Make output folder
outFolder="$inputsDir"/variantsCombo
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
#grep "$genotype" "$inputFileList" > tmpList_genotype.txt
#Add file type to end of each sample path
typeTag=_hap.g.vcf.gz
sed -e "s/$/$typeTag/" "$inputFileList" > tmpList.txt
#Add directory to beginning of each sample path
inDir="$inputsDir"/variantsCalled/
inDirTag=$(echo $inDir | sed "s/\//SLASH/g")
sed -i -e "s/^/$inDirTag/" tmpList.txt
#Add in slashes
sed -i "s/SLASH/\//g" tmpList.txt

#Add sample names to each path line
paste -d "$inputSampleList" tmpList.txt > tmpList_full.txt

#Output status mesasge
echo "Combining variants for the following input set of gvcf files: " > "$inputOutFile"
cat tmpList_full.txt >> "$inputOutFile"

#Combine gvcf files from multiple samples over the specified interval
for i in {1..496}; do
	echo "Combining variants for the following scaffold: $i"
	gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path "$outFolder" -L scaffold_"$i" --sample-name-map tmpList_full.txt --tmp-dir "$scratchPath" --reader-threads 4
	echo gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path "$outFolder" -L scaffold_"$i" --sample-name-map tmpList_full.txt --tmp-dir "$scratchPath" --reader-threads 4 >> "$inputOutFile"
done

#Clean up
rm tmpList*.txt
