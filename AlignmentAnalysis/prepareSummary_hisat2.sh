#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash prepareSummary_overall_hisat2.sh jobOutput alignmentMethod
#Usage Ex: bash prepareSummary_overall_hisat2.sh alignment_hisat2_jobOutput.o522510 hisat2
#Determine if the folder name was input in the correct format
if [[ $1 == *\/* ]] || [[ $1 == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Move to directory with output alignment folders
cd ../..
#Prepare input and output file names
inputStats="$1"
outputStats=alignmentStats_"$2"
#Retrieve sample order
grep "is being aligned" "$inputStats" > "$outputStats"_sampleOrder.txt
sed -i "s/.is being aligned\.\.\.//g" "$outputStats"_sampleOrder.txt
sed -i "s/Sample.//g" "$outputStats"_sampleOrder.txt
#Retrieve percent paired
#grep "were paired" "$inputStats" > "$outputStats"_paired.txt
#sed -i "s/; of these://g" "$outputStats"_paired.txt
#tr -s " " < "$outputStats"_paired.txt > tmp.txt"
#sed "s/^.//g" tmp.txt" > "$outputStats"_paired.txt
#Retrieve concordant alignment percents
grep "aligned concordantly" "$inputStats" > "tmp.txt"
grep -v "pairs aligned concordantly" "tmp.txt" > "$outputStats"_concordant.txt
tr -s " " < "$outputStats"_concordant.txt > "tmp.txt"
sed "s/^.//g" "tmp.txt" > "$outputStats"_concordant.txt
sed -i 'N;N;s/\n/,/g' "$outputStats"_concordant.txt
#Retrieve overall alignment percent
grep "overall alignment rate" "$inputStats" > "$outputStats"_overall.txt
#Combine all stats by sample
paste -d "," "$outputStats"_sampleOrder.txt "$outputStats"_concordant.txt "$outputStats"_overall.txt > "$outputStats"_combined.txt
#Clean up
rm "tmp.txt"
rm "$outputStats"_sampleOrder.txt
rm "$outputStats"_concordant.txt
rm "$outputStats"_overall.txt