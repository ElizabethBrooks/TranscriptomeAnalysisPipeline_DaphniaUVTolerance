#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_hisat2_sample.sh jobOutput alignmentMethod
#Usage Ex: bash alignmentSummary_hisat2_sample.sh alignment_hisat2_jobOutput.o522510 hisat2
#Move to directory with output alignment folders
cd ../..
#Prepare input and output file names
inputStats="$1"/alignedSummary.txt
outputStats=alignmentSummarized_"$2"
#Retrieve run number for input alignment folder
runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
#Retrieve sample order
grep "is being aligned" "$inputStats" > "$outputStats"_sampleOrder_"$runNum".txt
sed -i "s/.is being aligned\.\.\.//g" "$outputStats"_sampleOrder_"$runNum".txt
sed -i "s/Sample.//g" "$outputStats"_sampleOrder_"$runNum".txt
#Retrieve percent paired
#grep "were paired" "$inputStats" > "$outputStats"_paired_"$runNum".txt
#sed -i "s/; of these://g" "$outputStats"_paired_"$runNum".txt
#tr -s " " < "$outputStats"_paired_"$runNum".txt > tmp_"$runNum".txt"
#sed "s/^.//g" tmp_"$runNum".txt" > "$outputStats"_paired_"$runNum".txt
#Retrieve concordant alignment percents
grep "aligned concordantly" "$inputStats" > "tmp_"$runNum".txt"
grep -v "pairs aligned concordantly" "tmp_"$runNum".txt" > "$outputStats"_concordant_"$runNum".txt
tr -s " " < "$outputStats"_concordant_"$runNum".txt > "tmp_"$runNum".txt"
sed "s/^.//g" "tmp_"$runNum".txt" > "$outputStats"_concordant_"$runNum".txt
sed -i 'N;N;s/\n/,/g' "$outputStats"_concordant_"$runNum".txt
#Retrieve overall alignment percent
grep "overall alignment rate" "$inputStats" > "$outputStats"_overall_"$runNum".txt
#Combine all stats by sample
paste -d "," "$outputStats"_sampleOrder_"$runNum".txt "$outputStats"_concordant_"$runNum".txt "$outputStats"_overall_"$runNum".txt > "$outputStats"_combined_"$runNum".csv
#Clean up
rm "tmp_"$runNum".txt"
rm "$outputStats"_sampleOrder_"$runNum".txt
rm "$outputStats"_concordant_"$runNum".txt
rm "$outputStats"_overall_"$runNum".txt