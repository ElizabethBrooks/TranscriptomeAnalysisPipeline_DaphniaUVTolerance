#!/bin/bash
#Script to generate the longest CDS from a CDS fasta
#Usage: bash generateLongestCDS.sh inputFeatures inputReference genomeTag
#Usage Ex: bash generateLongestCDS.sh /Users/bamflappy/PfrenderLab/NCBI_dataset_Daphnia_pulex_Mar2022/data/KAP4_GCF_021134715.1/KAP4_cds.fa KAP4

#Retrieve input files
inputCDS=$1

#Set output file names
outDir=$(dirname $inputCDS)

#Create single line CDS file
tmpCDS="$outDir"/tmp_cds.fa
cat "$inputCDS" | sed -e ':a' -e 'N;$!ba' -e 's/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpCDS"

#Get list of CDS tags
tmpList="$outDir"/tmp_longest_cdsList.txt
cat "$inputCDS" | grep ">" | sed "s/ loc:.*//g" | sed "s/>//g" > "$tmpList"

#Loop over each CDS tag and retain longest CDS for each
outLongCDS=$outDir"/"$2"_longest_cds.fa"
outLongCDSList=$outDir"/"$2"_longest_cds_list.txt"
#Clear any previous data
[ -f $outLongCDS ] && rm $outLongCDS
echo "longestCDS" >> "$outLongCDSList"
#Output status messsage
echo "Writing longest CDS to file: $outLongCDS"
echo "Writing list of longest CDS to file: $outLongCDSList"
#Initialize cds length var
cLen=0
while IFS= read -r line; do
    #Check if current gene has multiple CDS ORF
    gTag=">$line"
    gLen=0
    gLong="$gTag"
    numCDS=$(grep -w "$gTag" "$tmpCDS" | wc -l)
    if [ $numCDS -gt 1 ]; then
        #Get the longest CDS
        for i in $(seq 1 $numCDS); do
            #Initialize variables
            cTag=$(grep -w "$gTag" "$tmpCDS" | head -$i | tail -1 | cut -d" " -f1)
            cStart=$(grep -w "$gTag" "$tmpCDS" | head -$i | tail -1 | grep "loc:*" | sed "s/segs:.*//g" | sed "s/^>.*)//g" | cut -d "-" -f1)
            cEnd=$(grep -w "$gTag" "$tmpCDS" | head -$i | tail -1 | grep "loc:*" | sed "s/segs:.*//g" | sed "s/^>.*)//g"| cut -d "-" -f2)
            cLen=$(($cEnd-$cStart))
            #Keep the longest CDS
            if [ $cLen -gt $gLen ]; then
                gLen=$cLen
                gLong="$cTag"
            fi
            cLen=0
        done
    fi
    #Output longest CDS
    grep -w "$gLong" "$tmpCDS" | sed 's/NEWLINE/\n/g' | sed "s/loc.*//g" >> "$outLongCDS"
    echo "$gLong" >> "$outLongCDSList"
done < "$tmpList"
echo "File with longest CDS has been generated!"

#Clean up
rm "$tmpCDS"
rm "$tmpList"