#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateKaKs_jobOutput

# script to generate MSAs for each gene in the reference set of peptide sequences
# usage: bash generateKaKs_musclePal2nalCodeml.sh
# usage ex: bash generateKaKs_musclePal2nalCodeml.sh

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# set software path
softwarePath=$(grep "pal2nal:" $baseDir"/InputData/softwarePaths.txt" | tr -d " " | sed "s/pal2nal://g")

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" $baseDir"/InputData/outputPaths.txt" | tr -d " " | sed "s/aligningGenome://g")
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
refPath=$(grep "genomeReference" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")
#refPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# retrieve file name of reference
refTag=$(basename $refPath)

# set reference pep and cds paths
inRefPep=$inputsPath"/"$refTag"_longest_pep.fa"
inRefNuc=$inputsPath"/"$refTag"_longest_cds.fa"

# set consensus pep and cds paths
inConPep=$inputsPath"/"$type"_consensus_longest_pep.fa"
inConNuc=$inputsPath"/"$type"_consensus_longest_cds.fa"

# set results file path
resultsFile=$outFolder"/kaksResults.csv"

# status message
echo "Begining analysis..."

# prepare reference multiline pep fasta to retrieve seqs
tmpRefPep=$outFolder"/Pulex_pep.tmp.fa"
cat $inRefPep | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpRefPep

# prepare input consensus multiline pep fasta
tmpConPep=$outFolder"/Olympics_pep.tmp.fa"
cat $inConPep | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpConPep

# prepare reference multiline cds fasta to retrieve seqs
tmpRefNuc=$outFolder"/Pulex_cds.tmp.fa"
cat $inRefNuc | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpRefNuc

# prepare input consensus multiline cds fasta
tmpConNuc=$outFolder"/Olympics_cds.tmp.fa"
cat $inConNuc | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpConNuc

# save ka ks values to final results file
resultsFile="$outFolder"/kaksResults.csv
echo "geneID  t  S  N  dNdS  dN  dS" > "$resultsFile"

# move to directory of inputs for running codeml
cd $outFolder

# loop over all genes in the reference
while IFS= read -r line; do
	# retrieve gene tag
	gTag=$(echo "$line" | sed 's/NEWLINE/\n/g' | grep ">" | cut -d " " -f 1 | sed 's/>//g')

	# set gene output paths
	gFile=$outFolder"/"$gTag"_Daphnia_pep.tmp.fa"
	gRefNuc=$outFolder"/"$gTag"_tmpRefNuc_cds.fa"
	gConNuc=$outFolder"/"$gTag"_tmpConNuc_cds.fa"
	outAln=$outFolder"/"$gTag"_Daphnia_aligned_pep.tmp.fa"

	# retrieve peptide sequences and convert back to multiline fasta format
	grep "^>$gTag" $tmpRefPep | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Pulex_$gTag/g" > $gFile
	grep "^>$gTag" $tmpConPep | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" >> $gFile

	# retrieve nucleotide sequences and convert back to multiline fasta format
	grep "^>$gTag" $tmpRefNuc | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Pulex_$gTag/g" > $gRefNuc
	grep "^>$gTag" $tmpConNuc | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Olympics_$gTag/g" > $gConNuc

	# prepare tree file
	echo "(>Pulex_$gTag, >Olympics_$gTag);" > $outFolder"/"$gTag".tree"

	# prepare control file template from pal2nal
	# seqtype = 1 for codon alignments
	# runmode = -2 performs ML estimation of dS and dN in pairwise comparisons
	cat $baseDir"/util/test.cnt" | sed "s/test\.codon/$gTag\.codon/g" | sed "s/test\.tree/$gTag\.tree/g" | sed "s/test\.codeml/$gTag\.codeml/g" > $outFolder"/"$gTag".cnt"

	# Status message
	echo "Generating MSA for $gTag..."

	# create MSA using muscle
	# https://stackoverflow.com/questions/70769809/muscle-command-line-wrapper
	muscle -in "$gFile" -out "$outAln"
	#muscle -align "$gFile" -output "$outAln"

	# status message
	echo "Generating codon alignment for $gTag..."

	# usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
	"$softwarePath"/pal2nal.pl "$outAln" "$gRefNuc" "$gConNuc" -output paml -nogap  >  $outFolder"/"$gTag".codon"
	#pal2nal.pl "$outAln" "$gRefNuc" "$gConNuc" -output paml -nogap  >  $outFolder"/"$gTag".codon"

	# status message
	echo "Generating ka ks values for $gTag..."

	# run codeml to retrieve ka ks values
	# you can find the output of codeml in the .codeml file
	# Ks, Ka values are very end of the output file
	codeml $gTag".cnt"

	# save ka ks values to final results file
	kaks=$(tail -1 $gTag".codeml")
	echo "$gTag  $kaks" >> $resultsFile

	# clean up
	rm $gFile
	rm $gRefNuc
	rm $gConNuc
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
rm $tmpRefPep
rm $tmpConPep
rm $tmpRefNuc
rm $tmpConNuc

# status message
echo "Analysis complete!"

