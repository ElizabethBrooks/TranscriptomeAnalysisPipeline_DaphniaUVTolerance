#!/bin/bash
#Script to generate the longest CDS from the PA42 v4.1 gff and reference fasta
#Usage: qsub generateLongestCDS_PA42_v4.1_cufflinks.sh

#Load necessary module
module load bio

#Retrieve input files
inputFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
inputRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Set output file names
outDir=$(dirname "$inputRef")
outCDS="$outDir"/PA42_v4.1_cds.fa

#Retrieve CDS
echo "Retrieving CDS..."
#Usage: gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] [-o <outfile.gff>] [-t <tname>] [-r #[[<strand>]<chr>:]<start>..<end> [-R]] [-CTVN‐ JMKQAFGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>] [-i <maxintron>]
gffread "$inputFeat" -g "$inputRef" -x "$outCDS" -W -F

#Create single line CDS file
tmpCDS="$outDir"/tmpPA42_v4.1_cds.fa
cat "$outCDS" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpCDS"

#Get list of CDS tags
tmpList="$outDir"/tmpPA42_v4.1_longest_cdsList.txt
colRefIn=$(grep "genePEPMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genePEPMap://g")
cat "$colRefIn" | cut -f1 > "$tmpList"

#Loop over each gene and retain longest CDS for each
outLongCDS="$outDir"/PA42_v4.1_longest_cds.fa
outLongCDSList="$outDir"/PA42_v4.1_longest_cds_list.txt
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
            cStart=$(grep -w "$gTag" "$tmpCDS" | head -$i | tail -1 | grep "loc:scaffold_*" | sed "s/segs:.*//g" | sed "s/^>.*)//g" | cut -d "-" -f1)
            cEnd=$(grep -w "$gTag" "$tmpCDS" | head -$i | tail -1 | grep "loc:scaffold_*" | sed "s/segs:.*//g" | sed "s/^>.*)//g"| cut -d "-" -f2)
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
    grep -w "$gLong" "$tmpCDS" | sed 's/NEWLINE/\n/g' | sed "s/loc:scaffold_.*//g" | sed "s/ gene=.*//g" >> "$outLongCDS"
    echo "$gLong" >> "$outLongCDSList"
done < "$tmpList"
echo "File with longest CDS has been generated!"

#Clean up
rm "$tmpCDS"
rm "$tmpList"