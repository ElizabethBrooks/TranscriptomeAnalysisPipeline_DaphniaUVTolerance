#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N annotation_trinotate_jobOutput
#Script to use trinotate to generate an annotation report
#Usage: qsub annotation_trinotate.sh transcriptomeFasta
#Usage Ex: qsub annotation_trinotate.sh trimmed_run1E05_assemblyTrinity
#Alternate usage Ex: qsub annotation_trinotate.sh PA42

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "trinotatePackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/trinotatePackage://g")
#Retrieve sqlite database path
sqliteDB=$(grep "trinotateSqliteDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/trinotateSqliteDB://g")
#Determine input query transcriptome for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/annotated_trinotate
elif [[ "$1" == PA42 ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
	#Set outputs absolute path
	inputsPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/annotated_trinotate
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Make output directory
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/annotated_trinotate_"$1"_summary.txt
#Load transcripts and coding regions
"$softsPath"/Trinotate "$sqliteDB" init --gene_trans_map "$inputsPath"/Trinity.fasta.gene_trans_map --transcript_fasta "$inputsPath"/Trinity.fasta --transdecoder_pep "$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
echo "$softsPath"/Trinotate "$sqliteDB" init --gene_trans_map "$inputsPath"/Trinity.fasta.gene_trans_map --transcript_fasta "$inputsPath"/Trinity.fasta --transdecoder_pep "$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.pep > "$inputOutFile"
#Loading BLAST homologies
"$softsPath"/Trinotate "$sqliteDB" LOAD_swissprot_blastp "$inputsPath"/searched_blastp_swissprot/blastp.outfmt6
 echo "$softsPath"/Trinotate "$sqliteDB" LOAD_swissprot_blastp "$inputsPath"/searched_blastp_swissprot/blastp.outfmt6 >> "$inputOutFile"
#Load Pfam domain entries
"$softsPath"/Trinotate "$sqliteDB" LOAD_pfam "$inputsPath"/searched_hmmscan/TrinotatePFAM.out
echo "$softsPath"/Trinotate "$sqliteDB" LOAD_pfam "$inputsPath"/searched_hmmscan/TrinotatePFAM.out >> "$inputOutFile"
#Optional load transmembrane domains
#"$softsPath"/Trinotate "$sqliteDB" LOAD_tmhmm tmhmm.out
#Optional load signal peptide predictions
#"$softsPath"/Trinotate "$sqliteDB" LOAD_signalp signalp.out
#Run Trinotate to generate an annotation report
"$softsPath"/Trinotate "$sqliteDB" report [opts] > trinotate_annotation_report.xls
echo "$softsPath"/Trinotate "$sqliteDB" report "[opts] >" trinotate_annotation_report.xls >> "$inputOutFile"
