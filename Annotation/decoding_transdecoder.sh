#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#$ -pe smp 1
#$ -q debug
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder
#Usage: qsub decoding_transdecoder.sh assembledTranscriptomeFolder
#Usage Ex: qsub decoding_transdecoder.sh trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98
#Usage Ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsHisat2_run2E05_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98
#Alternate usage Ex: qsub decoding_transdecoder.sh PA42_v4.1_cds
#Alternate usage Ex: qsub decoding_transdecoder.sh PA42_v4.1_transcripts
#Usage ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsHisat2_run3/variantCallingBcftools_filteredMapQ longest
#Usage ex: qsub decoding_transdecoder.sh genome

#Load necessary modules for ND CRC servers
module load bio
#module load bio/blast+
#module load bio/hmmer
#module load bio/cufflinks
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input query transcriptome for blastp
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve input assembly path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	#Set outputs absolute path
	outputsPath=$inputsPath/$1
	#Retrieve genome reference and features paths
	multiFASTA=$(echo "$outputsPath"/*.fasta)
	geneMap=$(echo "$outputsPath"/*.gene_trans_map)
	#Set output path
	outputFolder="$outputsPath"/"decoded_transdecoder"
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	#Set outputs absolute path
	outputsPath=$inputsPath/$1
	#Retrieve genome reference and features paths
	multiFASTA=$(echo "$outputsPath"/*.fasta)
	geneMap=$(echo "$outputsPath"/*.gene_trans_map)
	#Set output path
	outputFolder="$outputsPath"/"decoded_transdecoder"
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$inputsPath"
	#Retrieve genome reference and features paths
	geneMap=$(grep "geneCDSMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneCDSMap://g")
	#Set output path
	outputFolder="$outputsPath"/"decoded_transdecoder"
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$inputsPath"
	#Retrieve genome reference and features paths
	geneMap=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")
	#Set output path
	outputFolder="$outputsPath"/"decoded_transdecoder"
elif [[ "$1" == sorted* ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsPath="$inputsPath"/"$1"/transcripts_cufflinks.fa
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$inputsPath"
	#Retrieve genome reference and features paths
	geneMap="$outputsPath"/transcripts_cufflinks.fa.gene_trans_map
	#Set output path
	outputFolder="$outputsPath"/"decoded_transdecoder"
elif [[ "$1" == genome ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$outputsPath"/transcripts_cufflinks.fa
	#Retrieve genome reference and features paths
	geneMap="$outputsPath"/transcripts_cufflinks.fa.gene_trans_map
	#Set output path
	outputFolder="$outputsPath"/"decoded_transdecoder"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Determine if the input is clustered
if [[ "$1" == *clusteredNucleotide* ]]; then
	#Retrieve genome reference and features paths
	multiFASTA=$(echo "$outputsPath"/cdhitEst)
	inputsDir=$(dirname "$outputsPath")
	geneMap=$(echo "$inputsDir"/*.gene_trans_map)
elif [[ "$1" == *clusteredProtein* ]]; then
	#Retrieve genome reference and features paths
	multiFASTA=$(echo "$outputsPath"/cdhit)
	inputsDir=$(dirname "$outputsPath")
	geneMap=$(echo "$inputsDir"/*.gene_trans_map)
fi
#Make output folder
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile=$(echo "$1" | sed s'/\//./g')
inputOutFile="$outputFolder"/"$inputOutFile"_decoded_transdecoder_summary.txt

#Generate your best candidate open rading frame (ORF) predictions
#Generate candidate ORFs
echo "Beginning transdecoder open reading frame predictions..."
if [[ "$2" == "longest" || "$2" == "Longest" ]]; then
	TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap" --single_best_only
else
	TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap"
fi

#Output run commands to summary file
echo "TransDecoder.LongOrfs -t" "$multiFASTA" "--gene_trans_map" "$geneMap" > "$inputOutFile"
echo "Finished generating transdecoder open reading frame predictions!"
echo "Beginning transdecoder coding region selection..."
if [[ "$1" == PA42* || "$1" == sorted*  || "$1" == genome ]]; then
	TransDecoder.Predict -t "$multiFASTA" --no_refine_starts
	#Output run commands to summary file
	echo "TransDecoder.Predict -t" "$multiFASTA" "--no_refine_starts" >> "$inputOutFile"
else
	TransDecoder.Predict -t "$multiFASTA"
	#Output run commands to summary file
	echo "TransDecoder.Predict -t" "$multiFASTA" >> "$inputOutFile"
fi
echo "Finished transdecoder coding region selection!"
echo "Decoding complete!"
