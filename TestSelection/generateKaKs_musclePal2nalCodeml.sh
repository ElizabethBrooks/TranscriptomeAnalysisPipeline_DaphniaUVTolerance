#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateKaKs_jobOutput

# script to generate MSAs for each gene in the reference set of peptide sequences
# usage: bash generateKaKs_musclePal2nalCodeml.sh $tmpRefPep $tmpConPep $tmpRefNuc $tmpConNuc $resultsFile

# load necessary modules
module load bio

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# set software path
softwarePath=$(grep "pal2nal:" $baseDir"/InputData/softwarePaths.txt" | tr -d " " | sed "s/pal2nal://g")

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" $baseDir"/InputData/outputPaths.txt" | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# make outputs directory name
outFolder=$inputsPath"/selectionTests"

# retrieve subset tag
subsetTag=$1

# save ka ks values to final results file
resultsFile=$outFolder"/kaksResults_"$subsetTag".csv"
echo "geneID  t  S  N  dNdS  dN  dS" > "$resultsFile"

# prepare reference multiline pep fasta to retrieve seqs
tmpRefPep=$outFolder"/Pulex.pep.flt"$subsetTag

# prepare input consensus multiline pep fasta
tmpConPep=$outFolder"/Olympics.pep.flt"$subsetTag

# prepare reference multiline cds fasta to retrieve seqs
tmpRefNuc=$outFolder"/Pulex.cds.flt"$subsetTag

# prepare input consensus multiline cds fasta
tmpConNuc=$outFolder"/Olympics.cds.flt"$subsetTag

# status message
echo "Begining analysis..."

# move to directory of inputs for running codeml
cd $outFolder

# loop over all genes in the reference
while IFS= read -r line; do
	# retrieve gene tag
	gTag=$(echo "$line" | cut -f1 | sed 's/>//g')

	# set gene output paths
	gFile=$outFolder"/"$gTag"_Daphnia.pep.tmp.fa"
	gRefNuc=$outFolder"/"$gTag"_refNuc.cds.tmp.fa"
	gConNuc=$outFolder"/"$gTag"_conNuc.cds.tmp.fa"
	outAln=$outFolder"/"$gTag"_Daphnia_aligned.pep.tmp.fa"

	# retrieve peptide sequences and convert back to multiline fasta format
	echo ">Pulex_$gTag" > $gFile
	grep -w "^>$gTag" $tmpRefPep | cut -f2 | sed 's/\./X/g' >> $gFile
	echo ">Olympics_$gTag" >> $gFile
	grep -w "^>$gTag" $tmpConPep | cut -f2 | sed 's/\./X/g' >> $gFile

	# retrieve nucleotide sequences and convert back to multiline fasta format
	echo ">Pulex_$gTag" > $gRefNuc
	grep -w "^>$gTag" $tmpRefNuc | cut -f2 >> $gRefNuc
	echo ">Olympics_$gTag" > $gConNuc
	grep -w "^>$gTag" $tmpConNuc | cut -f2 >> $gConNuc

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
done < $tmpRefPep

# clean up
#rm $resultsFile
#rm $tmpRefPep
#rm $tmpConPep
#rm $tmpRefNuc
#rm $tmpConNuc

# move back to current directory
cd $currDir

# status message
echo "Analysis complete!"

