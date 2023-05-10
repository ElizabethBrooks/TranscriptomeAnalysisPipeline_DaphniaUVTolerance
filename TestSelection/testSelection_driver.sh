#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N testSelection_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub testSelection_driver.sh
# usage ex: qsub testSelection_driver.sh

# load necessary modules
#module load bio

# retrieve protein sequences
#bash retrieveFeatures_gffread.sh

# retreive protein coding sequence transcript names
#bash retrieve_proteinCoding.sh

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# set software path
softwarePath=$(grep "pal2nal:" $baseDir"/InputData/softwarePaths.txt" | tr -d " " | sed "s/pal2nal://g")

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" $baseDir"/InputData/softwarePaths.txt" | tr -d " " | sed "s/genomeFeatures://g")
#genomeFeatures="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff"

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" $baseDir"/InputData/outputPaths.txt" | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics"

# retrieve genome reference absolute path for alignment
refPath=$(grep "genomeReference" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")
#refPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/selectionTests"
mkdir $outFolder
#Check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outFolder directory already exsists... please remove before proceeding."
#	exit 1
#fi

# set inputs folder
inputsPath=$inputsPath"/features_gffread"

# retrieve file name of reference
refTag=$(basename $refPath)

# set results file path
resultsFile=$outFolder"/kaksResults.csv"
echo "geneID  t  S  N  dNdS  dN  dS" > "$resultsFile"

# set reference multiline pep fasta to retrieve seqs
fltRefPep=$inputsPath"/Pulex.pep.flt.fa"

# set input consensus multiline pep fasta
fltConPep=$inputsPath"/Olympics.pep.flt.fa"

# set reference multiline cds fasta to retrieve seqs
fltRefNuc=$inputsPath"/Pulex.cds.flt.fa"

# set input consensus multiline cds fasta
fltConNuc=$inputsPath"/Olympics.cds.flt.fa"

# split the gene sequence sets into segments
split -l 1000 $fltRefPep $outFolder"/Pulex.pep.flt.fa."
split -l 1000 $fltConPep $outFolder"/Olympics.pep.flt.fa."
split -l 1000 $fltRefNuc $outFolder"/Pulex.cds.flt.fa."
split -l 1000 $fltConNuc $outFolder"/Olympics.cds.flt.fa."

# loop over each segment
for i in $outFolder"/Pulex.pep.flt.fa."*; do
	# retrieve subset tag
	subsetTag=$(basename $i | sed 's/Pulex\.pep\.flt\.fa\.//g')
	# output status message
	echo "Starting analysis for subset $subsetTag ..."
	# generate Ka and Ks values for protein sequences
	qsub generateKaKs_musclePal2nalCodeml.sh $subsetTag
done

# wait
# https://stackoverflow.com/questions/11525214/wait-for-set-of-qsub-jobs-to-complete


# merge each of the ka ks results files


# fix formatting of the results file
#finalResults="$outFolder"/Pulex_Olympics_kaksResults.csv
#cat $resultsFile | sed "s/  /,/g" | sed "s/=,/=/g" | sed "s/ //g" | sed "s/dN\/dS=//g" | sed "s/dN=//g" | sed "s/dS=//g" | sed "s/t=//g" | sed "s/,S=//g" | sed "s/N=//g" > "$finalResults"

# clean up
#rm $resultsFile
#rm $geneList
#rm $transList

# status message
echo "Analysis complete!"
