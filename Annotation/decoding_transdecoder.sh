#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#$ -pe smp 12
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder
#Usage: qsub decoding_transdecoder.sh assembledTranscriptomeFolder blastDB optionalLongest
#Usage Ex: qsub decoding_transdecoder.sh trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98 ncbi
#Usage Ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsHisat2_run2E05_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 ncbi
#Alternate usage Ex: qsub decoding_transdecoder.sh PA42_v4.1_cds ncbi
#Alternate usage Ex: qsub decoding_transdecoder.sh PA42_v4.1_transcripts ncbi
#Usage ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsHisat2_run3/variantCallingBcftools_filteredMapQ ncbi longest
#Usage ex: qsub decoding_transdecoder.sh genome ncbi longest 

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
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	#Set outputs absolute path
	outputsPath=$inputsPath/$1
	#Retrieve genome reference and features paths
	multiFASTA=$(echo "$outputsPath"/*.fasta)
	geneMap=$(echo "$outputsPath"/*.gene_trans_map)
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$inputsPath"
	#Retrieve genome reference and features paths
	geneMap=$(grep "geneCDSMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneCDSMap://g")
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$inputsPath"
	#Retrieve genome reference and features paths
	geneMap=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")
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
elif [[ "$1" == genome ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	#Retrieve genome reference and features paths
	multiFASTA="$outputsPath"/transcripts_cufflinks.fa
	#Retrieve genome reference and features paths
	geneMap="$outputsPath"/transcripts_cufflinks.fa.gene_trans_map
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

#Determine input database for blastp
if [[ "$2" == "ncbi" ]]; then
	#Set slected database to ncbi
	blastpPath=$(grep "ncbiDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$2" == "uniprot" ]]; then
	#Set slected database to uniprot
	blastpPath=$(grep "uniprotDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/uniprotDB://g")
else
	#Error message
	echo "Invalid database selection of $2 entered (ncbi or uniprot only)... exiting!"
	exit 1
fi

#Retrieve genome reference and features paths
pfamPath=$(grep "pfamDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/pfamDB://g")

#Set output path
outputFolder="$outputsPath"/decoded_transdecoder_"$2"
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

#Generate candidate ORFs
echo "Beginning transdecoder open reading frame predictions..."
TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap"
#Output run commands to summary file
echo "TransDecoder.LongOrfs -t" "$multiFASTA" "--gene_trans_map" "$geneMap" > "$inputOutFile"
echo "Finished generating transdecoder open reading frame predictions!"

#Use BlastP to search a protein database
tag=$(basename "$multiFASTA")
echo "Beginning blastp protein database search..."
blastp -query "$tag".transdecoder_dir/longest_orfs.pep -db "$blastpPath"  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 12 > blastp.outfmt6
echo "Finished blastp protein database search!"

#Search the peptides for protein domains using Pfam
echo "Beginning hammscan search of peptides for protein domains..."
hmmscan --cpu 12 --domtblout pfam.domtblout "$pfamPath" "$tag".transdecoder_dir/longest_orfs.pep
echo "Finished hammscan search of peptides for protein domains!"

#Generate your best candidate open rading frame (ORF) predictions
#Combine the Blast and Pfam search results into coding region selection
echo "Beginning transdecoder coding region selection..."
if [[ "$1" == PA42* || "$1" == sorted*  || "$1" == genome ]]; then
	if [[ "$3" == "longest" || "$3" == "Longest" ]]; then
		TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts --single_best_only
		#Output run commands to summary file
		echo "TransDecoder.Predict -t" "$multiFASTA" "--retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts --single_best_only" >> "$inputOutFile"
	else
		TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts
		#Output run commands to summary file
		echo "TransDecoder.Predict -t" "$multiFASTA" "--retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts" >> "$inputOutFile"
	fi
else
	if [[ "$3" == "longest" || "$3" == "Longest" ]]; then
		TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only
		#Output run commands to summary file
		echo "TransDecoder.Predict -t $multiFASTA --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only" >> "$inputOutFile"
	else
		TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
		#Output run commands to summary file
		echo "TransDecoder.Predict -t" "$multiFASTA" "--retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6" >> "$inputOutFile"
	fi
fi
echo "Finished transdecoder coding region selection!"
echo "Decoding complete!"
