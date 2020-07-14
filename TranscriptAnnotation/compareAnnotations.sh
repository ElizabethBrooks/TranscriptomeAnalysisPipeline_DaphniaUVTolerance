#!/bin/bash

#Remove header
#tail -n +2 "$1" > tmpDataIn.csv

#Loop over annotations
while IFS=, read -r f1 f2
do
	#Check annotation sets
	if grep -iF "$f1,$f2" $2; then
		echo "$f1,$f2" >> "$1"_"$2"_matching.csv
	else
		echo "$f1,$f2" >> "$1"_unique.csv
	fi
done < $1

#Clean up
#rm tmp*.csv