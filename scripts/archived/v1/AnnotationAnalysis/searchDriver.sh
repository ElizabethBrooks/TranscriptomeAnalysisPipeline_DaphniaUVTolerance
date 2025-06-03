#!/bin/bash
#Script to perform sequence searches using a selected program for an input transcript data set
#Usage: bash searchDriver.sh method PA42Target assembledFolder sampleList

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

# set outputs path
outputFolder=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")

# set outputs absolute folder name
searchTag="$1"
outputFolder=$outputFolder"/reciprocalSearched_blastp_"$searchTag

#Set merged summary file name and header
if [[ "$1" == RBH ]]; then
	#Set output file name
	inputTag=$(echo $3 | sed 's/\//_/g')
	dbTag=$(echo $2 | sed 's/\//_/g')
	outFile="$outputFolder"/RBHB/"$inputTag"_"$dbTag"_blastp_summary.txt
	#Add header to output summary file
	echo "query,db,queryHits,dbHits,bestHits" > "$outFile"
elif [[ "$1" == consensus* ]]; then
	#Set output file names
	inputTag=$(echo $3 | sed 's/\//_/g')
	dbTag=$(echo $2 | sed 's/\//_/g')
	outFile="$outputFolder"/RBHB/"$inputTag"_"$dbTag"_blastp_consensusSummary.txt
	#Add header to output summary file
	echo "query,db,consensus,queryRBH,dbRBH,consensusRBH" > "$outFile"
else
	echo "Invalid analysis method entered... exiting!"
	exit 1
fi

#Initialize variables
counter=0
inputOutFile="$outputFolder"/RBHB/"$inputTag"_"$dbTag"_inputsSummary.txt
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
			qsub consensusRBH.sh "$inputFolder" "$i" "$2" "$outFile"
			#Write inputs to summary file
			echo "qsub consensusRBH.sh "$inputFolder" "$i" "$2" "$outFile >> $inputOutFile
		fi
	fi
	counter=$(($counter+1))
done

#Check if plotting was selected
if [ "$1" == plot ]; then
	#Load necessary modules
	module load R
	#Retrieve output file name
	outFile="$outputFolder"/reciprocalSearched_blastp/"$inputTag"_"$dbTag"_blastp_summary.txt
	#Usage: Rscript blastpStats_barPlots.r title blastpSummaryFile
	echo "Plotting blastp results for $3 and $2..."
	Rscript blastpStats_barPlots.r "$3" "$outFile"
	echo "Blastp results for $3 and $2 have been plotted!"
	#Write inputs to summary file
	echo "Rscript blastpStats_barPlots.r "$3" "$outFile >> $inputOutFile
fi
