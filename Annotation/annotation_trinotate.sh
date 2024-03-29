#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N annotation_trinotate_jobOutput
#Script to use trinotate to generate an annotation report
#Usage: qsub annotation_trinotate.sh transcriptomeFasta
#Usage Ex: qsub annotation_trinotate.sh trimmed_run1E05_assemblyTrinity
#Alternate usage Ex: qsub annotation_trinotate.sh PA42_v4.1_cds
#Alternate usage Ex: qsub annotation_trinotate.sh PA42_v3.0_transcripts
#Usage ex: qsub annotation_trinotate.sh sortedCoordinate_samtoolsHisat2_run3 filteredMapQ

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "trinotatePackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/trinotatePackage://g")
#Retrieve sqlite database path
sqliteDB=$(grep "trinotateSqlite:" ../InputData/databasePaths.txt | tr -d " " | sed "s/trinotateSqlite://g")
#Determine input query transcriptome for blastp
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsPath="$assemblyPath"/"$1"
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/annotated_trinotate
	#Set input paths
	geneTransMap="$inputsPath"/Trinity.fasta.gene_trans_map
	transcriptFasta="$inputsPath"/Trinity.fasta
	transdecoderPep="$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsPath="$assemblyPath"/"$1"
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/annotated_trinotate
	#Set input paths
	geneTransMap="$inputsPath"/Trinity.fasta.gene_trans_map
	transcriptFasta="$inputsPath"/Trinity.fasta
	transdecoderPep="$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	#Set outputs absolute path
	inputsPath=$(dirname "$inputsPath")
	outputFolder="$inputsPath"/annotated_trinotate
	#Set input paths
	geneTransMap=$(grep "geneCDSMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneCDSMap://g")
	transcriptFasta=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	transdecoderPep="$inputsPath"/decoded_transdecoder/*.transdecoder.pep
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	#Set outputs absolute path
	inputsPath=$(dirname "$inputsPath")
	outputFolder="$inputsPath"/annotated_trinotate
	#Set input paths
	geneTransMap=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")
	transcriptFasta=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	transdecoderPep="$inputsPath"/decoded_transdecoder/*.transdecoder.pep
elif [[ "$1" == sorted* ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsPath="$inputsPath"/"$1"/variantCallingBcftools_"$2"
	#Set outputs absolute path
	outputFolder="$inputsPath"/annotated_trinotate
	#Set input paths
	geneTransMap=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")
	transcriptFasta="$inputsPath"/transcripts_cufflinks.fa
	transdecoderPep="$inputsPath"/decoded_transdecoder/*.transdecoder.pep
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Set input paths
trinotateDB="$outputFolder"/Trinotate.sqlite
swissprotBlastpDB="$inputsPath"/searched_blastp_swissprot/blastp.outfmt6
pfamDB="$inputsPath"/searched_hmmscan/TrinotatePFAM.out
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
#Copy a new sqlite DB for the current transcript set
cp "$sqliteDB" "$outputFolder"
#Load transcripts and coding regions
"$softsPath"/Trinotate "$trinotateDB" init --gene_trans_map "$geneTransMap" --transcript_fasta "$transcriptFasta" --transdecoder_pep "$transdecoderPep"
echo "$softsPath"/Trinotate "$trinotateDB" init --gene_trans_map "$geneTransMap" --transcript_fasta "$transcriptFasta" --transdecoder_pep "$transdecoderPep" > "$inputOutFile"
#Load BLAST homologies
"$softsPath"/Trinotate "$trinotateDB" LOAD_swissprot_blastp "$swissprotBlastpDB"
echo "$softsPath"/Trinotate "$trinotateDB" LOAD_swissprot_blastp "$swissprotBlastpDB" >> "$inputOutFile"
#Load Pfam domain entries
"$softsPath"/Trinotate "$trinotateDB" LOAD_pfam "$pfamDB"
echo "$softsPath"/Trinotate "$trinotateDB" LOAD_pfam "$pfamDB" >> "$inputOutFile"
#Optional load transmembrane domains
#"$softsPath"/Trinotate "$trinotateDB" LOAD_tmhmm tmhmm.out
#Optional load signal peptide predictions
#"$softsPath"/Trinotate "$trinotateDB" LOAD_signalp signalp.out
#Run Trinotate to generate an annotation report
"$softsPath"/Trinotate "$trinotateDB" report > trinotate_annotation_report.xls
echo "$softsPath"/Trinotate "$trinotateDB" report ">" trinotate_annotation_report.xls >> "$inputOutFile"
