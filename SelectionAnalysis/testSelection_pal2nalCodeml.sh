#!/bin/bash
#Usae: bash testSelection_pal2nal.sh genotypeTag sampleSet variantCallingDir
#Usae ex: bash testSelection_pal2nal.sh dp_gene9-mRNA-1.p1 sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ

#Load necessary module
module load bio

#Set software path
softwarePath=$(grep "pal2nal:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/pal2nal://g")

#Set inputs path
inPath=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
inPath="$inPath"/daphniaMSA_PA42_v4.1_pep

#Set outputs path
outPath=$(grep "kaks:" ../InputData/outputPaths.txt | tr -d " " | sed "s/kaks://g")
outPath="$outPath"/daphniaKaks_PA42_v4.1

#Set genotype tag
gTag="$1"

#Retrieve input data
inAln="$inPath"/"$gTag"_pep_allDaphnia_aligned.fasta

#Retrieve input consensus cds
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="$inputsPath"/"$2"/"$3"
inConNuc="$inputsPath"/decoded_transdecoder/transcripts_cufflinks.fa.transdecoder.cds

#Retrieve input reference cds
refPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath=$(dirname $refPath)
inRefNuc="$refPath"/decoded_transdecoder/transcripts_cufflinks.fa.transdecoder.cds

#prepare consensus data files
tmpConNuc="$outPath"/"$gTag"_tmpConNuc.fa.cds
tmpConSeq="$outPath"/"$gTag"_tmpConSeq.fa.cds
cat "$inConNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpConNuc"
grep "^>$gTag" "$tmpConNuc" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" > "$tmpConSeq"
rm "$tmpConNuc"
#echo "Test consensus sequence: "
#cat "$tmpConSeq"

#prepare reference data files
tmpRefNuc="$outPath"/"$gTag"_tmpRefNuc.fa.cds
tmpRefSeq="$outPath"/"$gTag"_tmpRefSeq.fa.cds
cat "$inRefNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRefNuc"
grep "^>$gTag" "$tmpRefNuc" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$tmpRefSeq"
rm "$tmpRefNuc"
#echo "Test reference sequence: "
#cat "$tmpRefSeq"

#Prepare tree file
echo "(>PA42_v4.1_$gTag, >Olympics_$gTag);" > "$outPath"/"$gTag".tree
#cat "$outPath"/"$gTag".tree

#Prepare control file
cp "$softwarePath"/for_paml/test.cnt "$outPath"/"$gTag".cnt
sed -i "s/test\.codon/$gTag\.codon/g" "$outPath"/"$gTag".cnt
sed -i "s/test\.tree/$gTag\.tree/g" "$outPath"/"$gTag".cnt
sed -i "s/test\.codeml/$gTag\.codeml/g" "$outPath"/"$gTag".cnt
#cat "$outPath"/"$gTag".cnt

#Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
echo "Generating codon alignment for $gTag..."
"$softwarePath"/pal2nal.pl "$inAln" "$tmpConSeq" "$tmpRefSeq" -output paml  -nogap  >  "$outPath"/"$gTag".codon
#cat "$outPath"/"$gTag".codon

#Move to directory of inputs for codeml
cd "$outPath"

#Run codeml to retrieve ka ks values
#You can find the output of codeml in the .codeml file
#Ks, Ka values are very end of the output file
echo "Generating ka ks values for $gTag..."
codeml  "$gTag".cnt
#cat "$gTag".codeml
echo "Values of ka and ks have been generated!"

#Clean up
rm "$inAln"
rm "$tmpConSeq"
rm "$tmpRefSeq"
rm "$outPath"/"$gTag".codon
rm "$outPath"/"$gTag".tree
rm "$outPath"/"$gTag".cnt
