#!/bin/bash
#Script to perform sequence searches using a selected program for an input transcript data set
#Usage: bash searchDriver.sh method PA42Target assembledFolder sampleList
#Usage Ex: bash searchDriver.sh RBH PA42_v3.0_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash searchDriver.sh RBH PA42_v3.0_proteins trimmed_run1/clusteredNucleotides_cdhit_0.98 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash searchDriver.sh consensus PA42_v3.0_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash searchDriver.sh plot PA42_v3.0_proteins trimmed_run1
#Usage Ex: bash searchDriver.sh RBH PA42_v3.0_proteins sortedCoordinate_samtoolsHisat2_run2 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash searchDriver.sh RBH PA42_v4.1_proteins sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash searchDriver.sh RBH PA42_v4.1_proteins dnaRepair/Dmel_Svetec_2016 Dmel
#Usage Ex: bash searchDriver.sh RBH PA42_v4.1_proteins dnaRepair/Tcast_Guo_2019 Dmel

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
	inputTag=$(echo $3 | sed 's/\//_/g')
	outFile="$outPath"/RBHB/"$inputTag"_"$2"_blastp_summary.txt
	#Add header to output summary file
	echo "query,db,queryHits,dbHits,bestHits" > "$outFile"
elif [[ "$1" == consensus ]]; then
	#Set output file names
	inputTag=$(echo $3 | sed 's/\//_/g')
	outFile="$outPath"/RBHB/"$inputTag"_"$2"_blastp_consensusSummary.txt
	outFileUnique="$outPath"/RBHB/"$inputTag"_"$2"_blastp_uniqueRBH.txt
	#Add header to output summary file
	echo "query,db,queryRBH,dbRBH,consensusRBH,queryUnique" > "$outFile"
	echo "query,db,queryHit,dbHit" > "$outFileUnique"
else
	echo "Invalid analysis method entered... exiting!"
	exit 1
fi
#Initialize variables
counter=0
inputOutFile="$outPath"/RBHB/"$inputTag"_"$2"_inputsSummary.txt
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Determine what type of data folder was input
	if [[ "$3" == trimmed*/clustered* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			assemblyDir=$(echo $3 | cut -d '/' -f1)
			clusteredDir=$(echo $3 | cut -d '/' -f2)
			inputFolder=$(echo "$assemblyDir""$i"_assemblyTrinity/"$clusteredDir")
		fi
	elif [[ "$3" == sorted*/clustered* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			assemblyDir=$(echo $3 | cut -d '/' -f1)
			clusteredDir=$(echo $3 | cut -d '/' -f2)
			genomeTag=$(echo $2 | sed 's/_c.*//g' | sed 's/_p.*//g' | sed 's/_t.*//g')
			inputFolder=$(echo "$assemblyDir""$i"_assembly"$genomeTag"Trinity/"$clusteredDir")
		fi
	elif [[ "$3" == trimmed* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			inputFolder=$(echo "$3""$i"_assemblyTrinity)
		fi
	elif [[ "$3" == sorted* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			genomeTag=$(echo $2 | sed 's/_c.*//g' | sed 's/_p.*//g' | sed 's/_t.*//g')
			inputFolder=$(echo "$3""$i"_assembly"$genomeTag"Trinity)
		fi
	else
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			inputFolder="$3"
		fi
	fi
	#Check input option
	if [[ "$1" == RBH ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: bash searchRBH_blastp.sh transcriptomeFastaFolder
			echo "Merging $i blastp results for $3 and $2..."
			qsub searchRBH.sh "$inputFolder" "$i" "$2" "$outFile"
			#Write inputs to summary file
			echo "qsub searchRBH.sh "$inputFolder" "$i" "$2" "$outFile >> $inputOutFile
		fi
	elif [[ "$1" == consensus ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: qsub consensusRBH_blastp.sh transcriptomeFastaFolder
			echo "Generating consensus RBH of $i blastp results for $3 and $2..."
			qsub consensusRBH.sh "$inputFolder" "$i" "$2" "$outFile" "$outFileUnique"
			#Write inputs to summary file
			echo "qsub consensusRBH.sh "$inputFolder" "$i" "$2" "$outFile" "$outFileUnique >> $inputOutFile
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
	#Write inputs to summary file
	echo "Rscript blastpStats_barPlots.r "$3" "$outFile >> $inputOutFile
fi
