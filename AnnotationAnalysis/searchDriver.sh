#!/bin/bash
#Script to perform sequence searches using a selected program for an input transcript data set
#Usage: bash searchDriver.sh method PA42Target assembledFolder sampleList
#Usage Ex: bash searchDriver.sh RBH PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_transcripts PA42_proteins
#Usage Ex: bash searchDriver.sh consensus PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_transcripts PA42_proteins
#Usage Ex: bash searchDriver.sh plot PA42_proteins trimmed_run1

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#set output summary file path
outPath=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")
#Set merged summary file name and header
if [[ "$1" == RBH ]]; then
	#Set output file name
	outFile="$outPath"/reciprocalSearched_blastp/"$3"_"$2"_blastp_summary.txt
	#Add header to output summary file
	echo "query,db,queryHits,dbHits,bestHits" > "$outFile"
elif [[ "$1" == consensus ]]; then
	#Set output file names
	outFile="$outPath"/reciprocalSearched_blastp/"$3"_"$2"_blastp_consensusSummary.txt
	outFileUnique="$outPath"/reciprocalSearched_blastp/"$3"_"$2"_blastp_uniqueRBH.txt
	#Add header to output summary file
	echo "query,db,queryRBH,dbRBH,consensusRBH,queryUnique" > "$outFile"
	echo "query,db,queryHit,dbHit" > "$outFileUnique"
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
	if [[ "$1" == RBH ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: bash searchRBH_blastp.sh transcriptomeFastaFolder
			echo "Merging $i blastp results for $3 and $2..."
			qsub searchRBH_blastp.sh "$inputFolder" "$i" "$2" "$outFile"
		fi
	elif [[ "$1" == consensus ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: qsub consensusRBH_blastp.sh transcriptomeFastaFolder
			echo "Generating consensus RBH of $i blastp results for $3 and $2..."
			qsub consensusRBH_blastp.sh "$inputFolder" "$i" "$2" "$outFile" "$outFileUnique"
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
