#!/bin/bash
bash testSelection_pal2nal.sh

#Set software path
softwarePath=$(grep "pal2nal:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/pal2nal://g")

#Set inputs path
inPath=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
inPath="$outPath"/daphniaMSA_PA42_v4.1_pep

#Set outputs path
outPath=$(grep "kaks:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
outPath="$outPath"/daphniaKaks_PA42_v4.1_pep

#Check if the folder already exists
mkdir "$outPath"
if [ $? -ne 0 ]; then
	echo "The $outPath directory already exsists... please remove before proceeding."
	exit 1
fi


#Retrieve test data
testAln="dp_gene9_mRNA_1_p1_pep_allDaphnia_aligned.fasta"
testConNuc="/afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/PA42_v4.1/sortedCoordinate_samtoolsHisat2_run3/variantCallingBcftools_filteredMapQ/decoded_transdecoder/transcripts_cufflinks.fa.transdecoder.cds"
testRefNuc="/afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/decoded_transdecoder/transcripts_cufflinks.fa.transdecoder.cds"

#prepare consensus test data files
gTag="p_gene9-mRNA-1.p1"
tmpConNuc="$outPath"/tmpConNuc.fa.cds
tmpConSeq="$outPath"/tmpConSeq.fa.cds
cat "$testConNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpConNuc"
grep "^>$gTag" "$tmpConNuc" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" > "$tmpConSeq"
cat "$tmpConSeq"

#prepare reference test data files
tmpRefNuc="$outPath"/tmpRefNuc.fa.cds
tmpRefSeq="$outPath"/tmpRefSeq.fa.cds
cat "$testRefNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > "$tmpRefNuc"
grep "^>$gTag" "$tmpRefNuc" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>PA42_v4.1_$gTag/g" > "$tmpRefSeq"
cat "$tmpRefSeq"


#Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
cd "$softwarePath"
pal2nal.pl "$inPath"/"$testAln" "$tmpConSeq" "$tmpRefSeq" -output paml  -nogap  >  for_paml/test.codon

#Move to directory of inputs for codeml
#cd for_paml
ls for_paml

#Run codeml to retrieve ka ks values
#You can find the output of codeml in "test.codeml"
#Ks, Ka values are very end of the output file
#codeml  test.cnt
