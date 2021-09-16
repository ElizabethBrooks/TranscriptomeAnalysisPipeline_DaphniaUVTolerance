#!/bin/bash
#Script to run pal2nal and generate ka ks values
#Usae: bash testSelection_musclePal2nalCodeml.sh genotypeTag sampleSet variantCallingDir
#Usae ex: bash testSelection_musclePal2nalCodeml.sh dp_gene1-mRNA1 sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ

#Load necessary module
#module load bio

#Set software path
softwarePath=$(grep "pal2nal:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/pal2nal://g")

#Set inputs path
inPath=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
inPath="$inPath"/daphniaMSA_PA42_v4.1_pep

#Set outputs path
outPath=$(grep "kaks:" ../InputData/outputPaths.txt | tr -d " " | sed "s/kaks://g")
outPath="$outPath"/daphniaKaks_PA42_v4.1
failFile="$outPath"/failedGeneTags.txt

#Set genotype tag
gTag="$1"

#Prepare pep and alignment file names
gFile="$outPath"/tmp_pep_allDaphnia_"$gTag".fasta
gFileCleaned="$outPath"/tmpCleaned_pep_allDaphnia_"$gTag".fasta
outAln="$inPath"/"$gTag"_pep_allDaphnia_alignedOut.fasta
inAln="$inPath"/"$gTag"_pep_allDaphnia_aligned.fasta

#Retrieve input consensus data
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="$inputsPath"/"$2"/"$3"
inConNuc="$inputsPath"/Olympics_longest_cds.fa
#type=$(echo "$3" | cut -d"_" -f2)
#conFeat="$inputsPath"/"$type"_consensusFeatures.gff
#conSens=$(cat "$conFeat" | grep -w "$gTag" | grep "CDS" | cut -f7 | head -1 | sed "s/+/1/g" | sed "s/-/0/g")

#Retrieve input reference data
refPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath=$(dirname $refPath)
inRefNuc="$refPath"/PA42_v4.1_longest_cds.fa
#refFeat=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#refSens=$(cat "$refFeat" | grep -w "$gTag" | grep "CDS" | cut -f7 | head -1 | sed "s/+/1/g" | sed "s/-/0/g")

#Move to directory with translation script
cd ../Data

#Prepare single line reference data file
tmpRefNuc="$outPath"/"$gTag"_tmpRefNuc.fa.cds
#Retrieve reference CDS
echo ">PA42_v4.1_$gTag" > "$tmpRefNuc"
cat "$inRefNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' | grep -w "^>$gTag" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*//g" | tr 'a-z' 'A-Z' >> "$tmpRefNuc"
#Translate reference CDS to pep
echo ">PA42_v4.1_$gTag" > "$gFile"
#Rscript translateCDS_longestORF_seqinr.r "$tmpRefNuc" "$refSens" >> "$gFile"
Rscript translateCDS_longestORF_seqinr.r "$tmpRefNuc" >> "$gFile"
echo "" >> "$gFile"

#Prepare single line consensus data file
tmpConNuc="$outPath"/"$gTag"_tmpConNuc.fa.cds
#Retrieve consensus CDS
echo ">Olympics_$gTag" > "$tmpConNuc"
cat "$inConNuc" | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' | grep -w "^>$gTag" | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*//g" | tr 'a-z' 'A-Z' >> "$tmpConNuc"
#Translate consensus CDS to pep
echo ">Olympics_$gTag" >> "$gFile"
#Rscript translateCDS_longestORF_seqinr.r "$tmpConNuc" "$conSens" >> "$gFile"
Rscript translateCDS_longestORF_seqinr.r "$tmpConNuc" >> "$gFile"
echo "" >> "$gFile"

if cat "$gFile" | grep -q "NoValidORFFound" ; then
	#No valid ORF found for a sequence
	echo "$gTag" >> "$failFile"
else
	#Replace stop codon * with X wildcard
	sed 's/\*/X/g'  "$gFile" > "$gFileCleaned"
	rm "$gFile"

	#Output Status message
	echo "Generating MSA for $gTag..."

	#Create pep MSA
	muscle -in "$gFileCleaned" -out "$outAln"
	rm "$gFileCleaned"

	#Replace X wildcard with stop codon *
	sed 's/X/\*/g'  "$outAln" > "$inAln"
	rm "$outAln"

	#Output status message
	echo "MSA created for $gTag: $inAln"

	#Prepare tree file
	echo "(>PA42_v4.1_$gTag, >Olympics_$gTag);" > "$outPath"/"$gTag".tree

	#Prepare control file
	cp "$softwarePath"/for_paml/test.cnt "$outPath"/"$gTag".cnt
	sed -i "s/test\.codon/$gTag\.codon/g" "$outPath"/"$gTag".cnt
	sed -i "s/test\.tree/$gTag\.tree/g" "$outPath"/"$gTag".cnt
	sed -i "s/test\.codeml/$gTag\.codeml/g" "$outPath"/"$gTag".cnt

	#Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
	echo "Generating codon alignment for $gTag..."
	"$softwarePath"/pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  "$outPath"/"$gTag".codon

	#Move to directory of inputs for codeml
	cd "$outPath"

	#Run codeml to retrieve ka ks values
	#You can find the output of codeml in the .codeml file
	#Ks, Ka values are very end of the output file
	echo "Generating ka ks values for $gTag..."
	codeml  "$gTag".cnt
	echo "Values of ka and ks have been generated!"

	#Save ka ks values to final results file
	resultsFile=kaksResults.csv
	kaks=$(tail -1 "$gTag".codeml)
	echo "$gTag  $kaks" >> "$resultsFile"

	#Clean up
	[ -f "rst" ] && rm "rst"
	[ -f "rst1" ] && rm "rst1"
	[ -f "rub" ] && rm "rub"
	rm "$gTag".codon
	rm "$gTag".tree
	rm "$gTag".cnt
fi

#Clean up
rm "$inAln"
rm "$tmpConNuc"
rm "$tmpRefNuc"

