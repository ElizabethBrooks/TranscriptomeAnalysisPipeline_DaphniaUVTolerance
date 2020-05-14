#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoderPfamBlastp_jobOutput
#$ -pe smp 24
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder in conjunction with Pfam and Blastp results
#Usage: qsub decoding_transdecoderPfamBlastp.sh deNovoAssembledTranscriptomeFolder databaseSelection
#Usage Ex: qsub decoding_transdecoderPfamBlastp.sh trimmed_run1Sierra_assemblyTrinity ncbi
#Usage Ex: qsub decoding_transdecoderPfamBlastp.sh sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyGenomeTrinity ncbi
#Usage Ex: qsub decoding_transdecoderPfamBlastp.sh sortedCoordinate_samtoolsTophat2_run1Sierra_assemblyGenomeTrinity ncbi

#Load necessary modules for ND CRC servers
module load bio/transdecoder
module load bio/blast+
module load bio/hmmer
module load bio/cufflinks
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != *assembly* ]]; then
	echo "ERROR: The $1 folder of assembly files were not found... exiting"
	exit 1
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
#Retrieve input assembly path
inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
#Retrieve genome reference and features paths
pfamPath=$(grep "pfamDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/pfamDB://g")
multiFASTA=$(echo "$inputsPath"/"$1"/Trinity.fasta)
geneMap="$inputsPath"/"$1"/Trinity.fasta.gene_trans_map
#Retrieve outputs absolute path
outputsPath=$(grep "decoding:" ../InputData/outputPaths.txt | tr -d " " | sed "s/decoding://g")
outputFolder="$outputsPath"/decodedTransdecoder_"$1"
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_decodingTransdecoder_summary.txt
#Generate your best candidate open rading frame (ORF) predictions
echo "Beginning decoding..."
#Generate candidate ORFs
echo "Beginning transdecoder open reading frame predictions..."
TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap"
echo "Finished generating transdecoder open reading frame predictions!"
#Use BlastP to search a protein database
echo "Beginning blastp protein database search..."
blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep -db "$blastpPath"  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6
echo "Finished blastp protein database search!"
#Search the peptides for protein domains using Pfam
echo "Beginning hammscan search of peptides for protein domains..."
hmmscan --cpu 24 --domtblout pfam.domtblout "$pfamPath" Trinity.fasta.transdecoder_dir/longest_orfs.pep
echo "Finished hammscan search of peptides for protein domains!"
#Combine the Blast and Pfam search results into coding region selection
echo "Beginning transdecoder coding region selection..."
TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
echo "Finished transdecoder coding region selection!"
echo "Decoding complete!"
#Output run commands to summary file
echo "TransDecoder.LongOrfs -t" "$multiFASTA" "--gene_trans_map" "$geneMap" > "$inputOutFile"
echo "blastp -query" "Trinity.fasta.transdecoder_dir/longest_orfs.pep -db" "$blastpPath"  "-max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 >" "blastp.outfmt6" >> "$inputOutFile"
echo "hmmscan --cpu 24 --domtblout" "pfam.domtblout" "$pfamPath" "Trinity.fasta.transdecoder_dir/longest_orfs.pep" >> "$inputOutFile"
echo "TransDecoder.Predict -t" "$multiFASTA" "--retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6" >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"