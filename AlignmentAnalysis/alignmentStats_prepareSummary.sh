#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentStats_prepareSummary.sh jobOutput alignmentMethod
#Usage Ex: bash alignmentStats_prepareSummary.sh alignment_hisat2_jobOutput.o522510 hisat2
#Prepare input and output file names
inputStats="$1"
outputStats="alignmentStats_$2"
#Retrieve sample order
grep "is being aligned" "$inputStats.txt" > "$outputStats_sampleOrder.txt"
sed -i "s/.is being aligned\.\.\.//g" "$outputStats_sampleOrder.txt"
sed -i "s/Sample.//g" "$outputStats_sampleOrder.txt"
#Retrieve percent paired
#grep "were paired" "$inputStats.txt" > "$outputStats_paired.txt"
#sed -i "s/; of these://g" "$outputStats_paired.txt"
#tr -s " " < "$outputStats_paired.txt" > tmp.txt"
#sed "s/^.//g" tmp.txt" > "$outputStats_paired.txt"
#Retrieve concordant alignment percents
grep "aligned concordantly" "$inputStats.txt" > "tmp.txt"
grep -v "pairs aligned concordantly" "tmp.txt" > "$outputStats_concordant.txt"
tr -s " " < "$outputStats_concordant.txt" > "tmp.txt"
sed "s/^.//g" "tmp.txt" > "$outputStats_concordant.txt"
sed -i 'N;N;s/\n/,/g' "$outputStats_concordant.txt"
#Retrieve overall alignment percent
grep "overall alignment rate" "$inputStats.txt" > "$outputStats_overall.txt"
#Combine all stats by sample
paste -d "," "$outputStats_sampleOrder.txt" "$outputStats_concordant.txt" "$outputStats_overall.txt" > "$outputStats_combined.txt"
#Clean up
rm "tmp.txt"
rm "$outputStats_sampleOrder.txt"
rm "$outputStats_concordant.txt"
rm "$outputStats_overall.txt"