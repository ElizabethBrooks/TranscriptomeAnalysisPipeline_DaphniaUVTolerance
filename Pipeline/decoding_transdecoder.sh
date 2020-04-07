#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder
#Usage: qsub decoding_transdecoder.sh deNovoAssembledTranscriptomeFolder proteinDB
#Usage Ex: qsub decoding_transdecoder.sh trimmed_run1Sierra_assembly_Trinity uniprot-_dna+repair_+AND+organism_daphnia_isoformsIncluded.fasta
#Note that the genome version input is for output file naming purposes only
#Also note that uniprot databases may be downloaded from the UniprotKB search page,
# and Pfam databases may be downloaded with wget from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases
#ex: wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.full.gz

#Load necessary modules for ND CRC servers
module load bio/transdecoder
module load bio/blast+
module load bio/hmmer
#module load bio/cufflinks
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
if [[ "$1"  != trimmed*assembly_Trinity ]]; then
	echo "ERROR: The "$1" folder of trimmed fq.gz files were not found... exiting"
	exit 1
fi
#Retrieve genome reference and features paths
inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
uniprotPath=$(grep "uniprotDBs:" ../InputData/outputPaths.txt | tr -d " " | sed "s/uniprotDBs://g")
pfamPath=$(grep "pfamDB:" ../InputData/outputPaths.txt | tr -d " " | sed "s/pfamDB://g")
uniprotDB="$uniprotPath"/"$2"
multiFASTA=$(echo "$inputsPath"/"$1"/Trinity*.fasta)
geneMap="$inputsPath"/"$1"/Trinity.fasta.gene_trans_map
#Retrieve outputs absolute path
outputsPath=$(grep "decoding:" ../InputData/outputPaths.txt | tr -d " " | sed "s/decoding://g")
outputFolder="$outputsPath"/decoded_"$1"
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_decoding_summary.txt
#Generate your best candidate open rading frame (ORF) predictions
echo "Beginning decoding..."
#Generate candidate ORFs
#TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap"
#Use BlastP to search a protein database
blastp -query "$outputFolder"/Trinity.fasta.transdecoder_dir/longest_orfs.pep -db "$uniprotDB"  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > "$outputFolder"/blastp.outfmt6
#Search the peptides for protein domains using Pfam
#hmmscan --cpu 8 --domtblout "$outputFolder"/pfam.domtblout "$pfamPath" "$outputFolder"/Trinity.fasta.transdecoder_dir/longest_orfs.pep
#Combine the Blast and Pfam search results into coding region selection
#TransDecoder.Predict -t "$multiFASTA" --retain_pfam_hits "$outputFolder"/pfam.domtblout --retain_blastp_hits "$outputFolder"/blastp.outfmt6
echo "Decoding finished!"
#Clean up
rm "$outputFolder"/blastp.outfmt6
rm "$outputFolder"/pfam.domtblout
#Output run commands to summary file
echo "TransDecoder.LongOrfs -t" "$multiFASTA" "--gene_trans_map" "$geneMap" > "$inputOutFile"
echo "blastp -query" "$outputFolder"/"Trinity.fasta.transdecoder_dir/longest_orfs.pep -db" "$uniprotDB"  "-max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 >" "$outputFolder"/"blastp.outfmt6" >> "$inputOutFile"
echo "hmmscan --cpu 8 --domtblout" "$outputFolder"/"pfam.domtblout" "$pfamPath" "$outputFolder"/"Trinity.fasta.transdecoder_dir/longest_orfs.pep" >> "$inputOutFile"
echo "TransDecoder.Predict -t" "$multiFASTA" >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"