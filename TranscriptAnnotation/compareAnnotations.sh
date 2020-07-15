#!/bin/bash
#Script to compare annotation files

#Remove header
#tail -n +2 "$1" > tmpDataIn.csv

#Loop over annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -iF "$f1,$f2" $2; then #Matching
		echo "$f1,$f2" >> "$1"_"$2"_matching.csv
	else #Unique
		echo "$f1,$f2" >> "$1"_unique.csv
	fi
done < $1

#Output the number of annotations in each subset
wc -l "$1"_"$2"_matching.csv
wc -l "$1"_unique.csv

#Clean up
#rm tmp*.csv