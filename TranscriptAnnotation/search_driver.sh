#!/bin/bash
#Script to perform sequence searches using a selected program for an input transcript data set
#Usage: bash search_driver.sh method PA42Target assembledFolder sampleList
#Usage Ex: bash search_driver.sh blastp PA42_proteins trimmed_run1 swissprot Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh hmmscan PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh reciprocal PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_transcripts PA42_proteins
#Usage Ex: bash search_driver.sh merge PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_transcripts PA42_proteins
#Usage Ex: bash search_driver.sh merge PA42_transcripts trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_proteins
#Usage Ex: bash search_driver.sh merge PA42_cds trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_proteins PA42_transcripts
#Usage Ex: bash search_driver.sh consensus PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_transcripts
#Usage Ex: bash search_driver.sh plot PA42_proteins trimmed_run1

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#set output summary file path
outPath=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")
#Set merged summary file name and header
if [[ "$1" == merge ]]; then
	#Set output file name
	outFile="$outPath"/reciprocalSearched_blastp/"$3"_"$2"_blastp_summary.txt
	#Add header to output summary file
	echo "query,db,queryHits,dbHits,bestHits" > "$outFile"
fi
#Initialize variables
counter=0
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Determine what type of data folder was input
	if [[ "$3" == trimmed* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			inputFolder=$(echo "$3""$i"_assemblyTrinity)
		fi
	elif [[ "$3" == sorted* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			inputFolder=$(echo "$3""$i"_assemblyGenomeTrinity)
		fi
	else
		echo "ERROR: Input folder for analysis is not a valid option... exiting!"
		exit 1
	fi
	#Check input option
	if [[ "$1" == reciprocal ]]; then #Skip first three arguments 
		if [ $counter -ge 3 ]; then
			#Usage: qsub reciprocalSearch_blastp.sh transcriptomeFastaFolder
			echo "Searching $i transcriptome with blastp for $3 and $2..."
			qsub reciprocalSearch_blastp.sh "$inputFolder" "$2"
		fi
	elif [[ "$1" == merge ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: bash mergeSearches_blastp.sh transcriptomeFastaFolder
			echo "Merging $i blastp results for $3 and $2..."
			qsub mergeSearches_blastp.sh "$inputFolder" "$i" "$2" "$outFile"
		fi
	elif [[ "$1" == consensus ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: bash consensusRBH_blastp.sh transcriptomeFastaFolder
			echo "Generating consensus RBH of $i blastp results for $3 and $2..."
			qsub consensusRBH_blastp.sh "$inputFolder" "$i" "$2" "$outFile"
		fi
	elif [[ "$1" == hmmscan ]]; then #Skip first three arguments 
		if [ $counter -ge 3 ]; then
			#Usage: qsub search_hmmscan.sh transcriptomeFastaFolder
			echo "Searching $i transcriptome with hmmscan for $3 and $2..."
			qsub search_hmmscan.sh "$inputFolder" "$2"
		fi
	elif [[ "$1" == blastp ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: qsub search_blastp.sh transcriptomeFastaFolder proteinDB
			echo "Searching $i transcriptome with blastp using $2..."
			qsub search_blastp.sh "$inputFolder" "$2"
		fi
	fi
	counter=$(($counter+1))
done
#Check if plotting was selected
if [ "$1" == plot ]; then
	#Load necessary modules
	module load R
	#Retrieve output file name
	outFile="$outPath"/reciprocalSearched_blastp/"$3"_"$2"_blastp_summary.txt
	#Usage: Rscript blastpStats_barPlots.r title blastpSummaryFile
	echo "Plotting blastp results for $3 and $2..."
	Rscript blastpStats_barPlots.r "$3" "$outFile"
	echo "Blastp results for $3 and $2 have been plotted!"
fi
