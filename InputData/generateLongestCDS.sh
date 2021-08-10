#!/bin/bash
#Script to generate the longest CDS from the PA42 v4.1 gff and reference fasta
#Usage: bash generateLongestCDS.sh

#Load necessary module
module load bio

#Retrieve input files
inputFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
inputRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Set output file names
outDir=$(dirname "$inputRef")
outCDS="$outDir"/PA42_v4.1_CDS.fa

#Usage: gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] [-o <outfile.gff>] [-t <tname>] [-r #[[<strand>]<chr>:]<start>..<end> [-R]] [-CTVN‐ JMKQAFGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>] [-i <maxintron>]
gffread "$inputFeat" -g "$inputRef" -x "$outCDS" -W -F

#Create single line CDS file
tmpCDS="$outDir"/tmpPA42_v4.1_longestCDS.fa
cat "$outCDS" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpCDS"

#Get list of CDS tags
tmpList="$outDir"/tmpPA42_v4.1_longestCDSList.txt
colRefIn=$(grep "genePEPMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genePEPMap://g")
cat "$colRefIn" | cut -f1 > "$tmpList"

#Loop over each gene and retain longest CDS for each
outLongCDS="$outDir"/PA42_v4.1_longestCDS.fa
[ -f $outLongCDS ] && rm $outLongCDS
cLen=0
while IFS= read -r line; do
    #Check if current gene has multiple CDS ORF
    gTag=">$line"
    gLen=0
    numCDS=$(grep -w "$gTag" "$tmpCDS" | wc -l)
    if [ $numCDS -gt 1 ]; then
        #Get the length of each CDS
        for i in $(seq 1 $numCDS); do
            #Initialize variables
            cTag=$(grep -w "$gTag" "$tmpCDS" | head -$i | tail -1 | cut -d" " -f1)
            cStart=$(grep -w "$cTag" "$tmpCDS" | head -$i | tail -1 | cut -d" " -f3 | cut -d ")" -f2 | cut -d "-" -f2)
            cEnd=$(grep -w "$cTag" "$tmpCDS" | head -$i | tail -1 | cut -d" " -f3 | cut -d ")" -f2 | cut -d "-" -f2)
            cLen=$(($gEnd-$gStart))
            #Keep the longest CDS
            if [ $cLen -gt $gLen ]; then
                gLen=$cLen
                gTag="$cTag"
            fi
            cLen=0
        done
    fi
    #Output longest CDS
    grep -w "$gTag" "$tmpCDS" | sed 's/NEWLINE/\n/g' >> "$outLongCDS"
done < "$tmpList"

#Clean up
rm "$tmpCDS"
rm "$tmpList"