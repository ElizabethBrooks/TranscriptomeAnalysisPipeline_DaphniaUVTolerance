#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_tophat2_sample.sh alignedSampleFolder alignmentMethod runNum
#Usage Ex: bash alignmentSummary_tophat2_sample.sh alignment_tophat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_Y05_UV topaht2 run2
#Prepare input and output file names
runNum="$3"
inputStats="$1"align_summary.txt
outputStats=alignmentSummarized_"$2"
#Retrieve run number for input alignment folder
runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
#Retrieve mapped left read percents
grep "Mapped" "$inputStats" > "$outputStats"_mapped_"$runNum".txt
head -1 "$outputStats"_mapped_"$runNum".txt > "$outputStats"_mappedLeft_"$runNum".txt
cat "$outputStats"_mappedLeft_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_mappedLeft_"$runNum".txt
sed -i "s/(//g" "$outputStats"_mappedLeft_"$runNum".txt
#Retrieve mapped right read percents
tail -1 "$outputStats"_mapped_"$runNum".txt > "$outputStats"_mappedRight_"$runNum".txt
cat "$outputStats"_mappedRight_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_mappedRight_"$runNum".txt
sed -i "s/(//g" "$outputStats"_mappedRight_"$runNum".txt
#Retrieve overall percents
grep "overall" "$inputStats" > "$outputStats"_overall_"$runNum".txt
cat "$outputStats"_overall_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_overall_"$runNum".txt
#Retrieve concordant percent
grep "concordant" "$inputStats" > "$outputStats"_concordant_"$runNum".txt
cat "$outputStats"_concordant_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_concordant_"$runNum".txt
#Combine all stats by sample
paste -d "," "$outputStats"_mapped_"$runNum".txt "$outputStats"_overall_"$runNum".txt "$outputStats"_concordant_"$runNum".txt > "$outputStats"_combined_"$runNum".csv
#Clean up
rm "tmp_"$runNum".txt"
rm "$outputStats"_mapped_"$runNum".txt
rm "$outputStats"_overall_"$runNum".txt
rm "$outputStats"_concordant_"$runNum".txt