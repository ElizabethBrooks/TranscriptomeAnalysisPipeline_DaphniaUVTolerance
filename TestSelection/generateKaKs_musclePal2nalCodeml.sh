#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateKaKs_jobOutput

# script to generate MSAs for each gene in the reference set of peptide sequences
# usage: bash generateKaKs_musclePal2nalCodeml.sh
# usage ex: bash generateKaKs_musclePal2nalCodeml.sh

# set software path
softwarePath=$(grep "pal2nal:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/pal2nal://g")

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/selectionTests"
mkdir $outFolder
#Check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outFolder directory already exsists... please remove before proceeding."
#	exit 1
#fi

# set inputs folder
inputsPath=$inputsPath"/features_gffread"

# retrieve genome reference absolute path for alignment
refPath=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#refPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# retrieve file name of reference
refTag=$(basename $refPath)

# set reference pep path
refPath=$inputsPath"/"$refTag"_longest_pep.fa"

# set consensus pep path
inputsPath=$inputsPath"/"$type"_consensus_longest_pep.fa"

# set results file path
resultsFile=$outFolder"/kaksResults.csv"

# status message
echo "Generating MSA for protein sequences..."

# prepare reference multiline pep fasta to retrieve seqs
tmpRef=$outFolder"/tmpPulex_pep.fa"
cat $refPath | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpRef

# prepare input multiline pep fasta
tmpSample=$outFolder"/tmpOLYM_pep.fa"
cat $inputsPath | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpSample

# retrieve base directory path
baseDir=$(dirname $currDir)

# save ka ks values to final results file
resultsFile="$outFolder"/kaksResults.csv
echo "geneID  t  S  N  dNdS  dN  dS" > "$resultsFile"

# loop over all genes in the reference
while IFS= read -r line; do
	# retrieve gene tag
	gTag=$(echo $line | sed 's/NEWLINE/\n/g' | grep ">" | cut -d " " -f 1 | sed 's/>rna-//g')

	# set tmp output pep file
	gFile=$outFolder"/tmp_"$gTag"_Daphnia_pep.fa"

	# retrieve selected peptide sequences and convert back to multiline fasta format
	grep "^>$gTag" $tmpRef | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Pulex_$gTag/g" > $gFile
	
	# Status message
	echo "Generating MSA for $gTag..."

	# prepare multiline pep fasta to retrieve seqs
	grep "^>$gTag" $tmpSample | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>OLYM_$gTag/g" >> $gFile

	# set file paths
	outAln=$outFolder"/tmp_"$gTag"_Daphnia_aligned_pep.fa"
	inAln=$outFolder"/"$gTag"_Daphnia_pep.fa"

	# create MSA using muscle
	# https://stackoverflow.com/questions/70769809/muscle-command-line-wrapper
	muscle -in "$gFile" -out "$outAln"
	#muscle -align "$gFile" -output "$outAln"

	# replace X wildcard with stop codon *
	sed 's/X/\*/g'  "$outAln" > "$inAln"

	# clean up
	rm $gFile
	rm $outAln

	# prepare tree file
	echo "(>Pulex_$gTag, >Olympics_$gTag);" > $outFolder"/"$gTag".tree"

	# prepare control file template from pal2nal
	# seqtype = 1 for codon alignments
	# runmode = -2 performs ML estimation of dS and dN in pairwise comparisons
	cat $baseDir"/util/test.cnt" | sed "s/test\.codon/$gTag\.codon/g" | sed "s/test\.tree/$gTag\.tree/g" | sed "s/test\.codeml/$gTag\.codeml/g" > $outFolder"/"$gTag".cnt"

	# status message
	echo "Generating codon alignment for $gTag..."

	# usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
	"$softwarePath"/pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  $outFolder"/"$gTag".codon"
	#pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  $outFolder"/"$gTag".codon"

	# move to directory of inputs for codeml
	cd $outFolder

	# status message
	echo "Generating ka ks values for $gTag..."

	# run codeml to retrieve ka ks values
	# you can find the output of codeml in the .codeml file
	# Ks, Ka values are very end of the output file
	codeml $outFolder"/"$gTag".cnt"

	# save ka ks values to final results file
	kaks=$(tail -1 $gTag".codeml")
	echo "$gTag  $kaks" >> $resultsFile

	# clean up
	[ -f "rst" ] && rm "rst"
	[ -f "rst1" ] && rm "rst1"
	[ -f "rub" ] && rm "rub"
	rm $gTag".codon"
	rm $gTag".tree"
	rm $gTag".cnt"
done < $tmpRef

# fix formatting of the results file
finalResults="$outFolder"/Pulex_Olympics_kaksResults.csv
cat $resultsFile | sed "s/  /,/g" | sed "s/=,/=/g" | sed "s/ //g" | sed "s/dN\/dS=//g" | sed "s/dN=//g" | sed "s/dS=//g" | sed "s/t=//g" | sed "s/,S=//g" | sed "s/N=//g" > "$finalResults"

# clean up
rm $tmpSample
rm $tmpRef

# status message
echo "MSA generated for protein sequences!"

