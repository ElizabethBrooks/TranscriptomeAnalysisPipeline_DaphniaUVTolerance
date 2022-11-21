#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transcripts_transdecoder_jobOutput
#$ -pe smp 12
# script to predict coding regions from transcript sequences using Transdecoder
# usage: qsub decoding_transcripts_transdecoder.sh blastDB optionalLongest
# usage Ex: qsub decoding_transcripts_transdecoder.sh ncbi
# usage Ex: qsub decoding_transcripts_transdecoder.sh ncbi longest

#Load necessary module for ND CRC servers
module load bio

# retrieve genome tag
genomeTag=$(grep "genomeTag:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeTag://g")
# retrieve genome reference absolute path for querying
inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
# retrieve genome reference and features paths
multiFASTA="$inputsPath"
# retrieve genome reference and features paths
geneMap=$(grep "geneTransMap:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTransMap://g")

# set outputs absolute path
outputsPath=$(grep "translation:" ../InputData/outputPaths.txt | tr -d " " | sed "s/translation://g")

#Determine input database for blastp
if [[ "$1" == "ncbi" ]]; then
	#Set slected database to ncbi
	blastpPath=$(grep "ncbi:" ../InputData/databasePaths.txt | tr -d " " | sed "s/ncbi://g")
elif [[ "$1" == "uniprot" ]]; then
	#Set slected database to uniprot
	blastpPath=$(grep "uniprot:" ../InputData/databasePaths.txt | tr -d " " | sed "s/uniprot://g")
else
	#Error message
	echo "Invalid database selection of $2 entered (ncbi or uniprot only)... exiting!"
	exit 1
fi

#Retrieve genome reference and features paths
pfamPath=$(grep "pfam:" ../InputData/databasePaths.txt | tr -d " " | sed "s/pfam://g")

#Set output path
if [[ "$2" == "longest" || "$2" == "Longest" ]]; then
	outputFolder="$outputsPath"/transcriptsDecoded_transdecoder_longest_"$genomeTag"
else
	outputFolder="$outputsPath"/transcriptsDecoded_transdecoder_"$genomeTag"
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
if [[ "$2" == "longest" || "$2" == "Longest" ]]; then
	TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts --single_best_only
	#Output run commands to summary file
	echo "TransDecoder.Predict -t" "$multiFASTA" "--retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts --single_best_only" >> "$inputOutFile"
else
	TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts
	#Output run commands to summary file
	echo "TransDecoder.Predict -t" "$multiFASTA" "--retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --no_refine_starts" >> "$inputOutFile"
fi
echo "Finished transdecoder coding region selection!"
echo "Decoding complete!"
