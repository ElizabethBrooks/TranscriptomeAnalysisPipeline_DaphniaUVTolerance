#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash prepareSummary_sample_tophat2.sh alignedSampleFolder alignmentMethod
#Usage Ex: bash prepareSummary_sample_tophat2.sh alignment_topaht2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_Y05_UV topaht2
#Determine if the folder name was input in the correct format
if [[ $1 == *\/* ]] || [[ $1 == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Prepare input and output file names
inputStats="$1"/align_summary.txt
outputStats=alignmentStats_"$2"
#Retrieve mapped left read percents
grep "Mapped" "$inputStats" > "$outputStats"_mapped.txt
head -1 "$outputStats"_mapped.txt > "$outputStats"_mappedLeft.txt
cat "$outputStats"_mappedLeft.txt | tr " " "\n" > tmp.txt
grep "%" tmp.txt > "$outputStats"_mappedLeft.txt
sed -i "s/(//g" "$outputStats"_mappedLeft.txt
#Retrieve mapped right read percents
tail -1 "$outputStats"_mapped.txt > "$outputStats"_mappedRight.txt
cat "$outputStats"_mappedRight.txt | tr " " "\n" > tmp.txt
grep "%" tmp.txt > "$outputStats"_mappedRight.txt
sed -i "s/(//g" "$outputStats"_mappedRight.txt
#Retrieve overall percents
grep "overall" "$inputStats" > "$outputStats"_overall.txt
cat "$outputStats"_overall.txt | tr " " "\n" > tmp.txt
grep "%" tmp.txt > "$outputStats"_overall.txt
#Retrieve concordant percent
grep "concordant" "$inputStats" > "$outputStats"_concordant.txt
cat "$outputStats"_concordant.txt | tr " " "\n" > tmp.txt
grep "%" tmp.txt > "$outputStats"_concordant.txt
#Combine all stats by sample
paste -d "," "$outputStats"_mapped.txt "$outputStats"_overall.txt "$outputStats"_concordant.txt > "$outputStats"_combined.csv
#Clean up
rm "tmp.txt"
rm "$outputStats"_mapped.txt
rm "$outputStats"_overall.txt
rm "$outputStats"_concordant.txt