#!/bin/bash
#bash testSelection_pal2nal.sh

#Set software path
softwarePath=$(grep "pal2nal:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/pal2nal://g")

#Set inputs path
inPath=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
inPath="$outPath"/daphniaMSA_PA42_v4.1_pep

#Set outputs path
outPath=$(grep "kaks:" ../InputData/outputPaths.txt | tr -d " " | sed "s/kaks://g")
outPath="$outPath"/daphniaKaks_PA42_v4.1_pep

#Check if the folder already exists
mkdir "$outPath"
if [ $? -ne 0 ]; then
	echo "The $outPath directory already exsists... please remove before proceeding."
	exit 1
fi


#Retrieve input data
inAln="dp_gene9_mRNA_1_p1_pep_allDaphnia_aligned.fasta"
inConNuc="/afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/PA42_v4.1/sortedCoordinate_samtoolsHisat2_run3/variantCallingBcftools_filteredMapQ/decoded_transdecoder/transcripts_cufflinks.fa.transdecoder.cds"
inRefNuc="/afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/decoded_transdecoder/transcripts_cufflinks.fa.transdecoder.cds"

#Set genotype tag
gTag="dp_gene9-mRNA-1.p1"

#prepare consensus data files
tmpConNuc="$outPath"/tmpConNuc.fa.cds
tmpConSeq="$outPath"/tmpConSeq.fa.cds
cat "$inConNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpConNuc"
grep "^>$gTag" "$tmpConNuc" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" > "$tmpConSeq"
rm "$tmpConNuc"
echo "Test consensus sequence: "
cat "$tmpConSeq"

#prepare reference data files
tmpRefNuc="$outPath"/tmpRefNuc.fa.cds
tmpRefSeq="$outPath"/tmpRefSeq.fa.cds
cat "$inRefNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRefNuc"
grep "^>$gTag" "$tmpRefNuc" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$tmpRefSeq"
rm "$tmpRefNuc"
echo "Test reference sequence: "
cat "$tmpRefSeq"

#Prepare tree file
echo "(>PA42_v4.1_$gTag, >Olympics_$gTag);" > "$outPath"/"$gTag".tree
cat "$outPath"/"$gTag".tree

#Prepare control file
cp "$softwarePath"/for_paml/test.cnt "$outPath"/"$gTag".cnt
sed -i "s/test\.codon/$gTag\.codon/g" "$outPath"/"$gTag".cnt
sed -i "s/test\.tree/$gTag\.tree/g" "$outPath"/"$gTag".cnt
sed -i "s/test\.codeml/$gTag\.codeml/g" "$outPath"/"$gTag".cnt
cat "$outPath"/"$gTag".cnt

#Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
"$softwarePath"/pal2nal.pl "$inPath"/"$inAln" "$tmpConSeq" "$tmpRefSeq" -output paml  -nogap  >  "$outPath"/"$gTag".codon
cat "$outPath"/"$gTag".codon

#Clean up
rm "$tmpConSeq"
rm "$tmpRefSeq"

#Move to directory of inputs for codeml
cd "$outPath"

#Run codeml to retrieve ka ks values
#You can find the output of codeml in the .codeml file
#Ks, Ka values are very end of the output file
"$softwarePath"/codeml  "$gTag".cnt
cat "$gTag".codeml
