#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateKaKs_jobOutput

# script to generate MSAs for each gene in the reference set of peptide sequences
# then create the codon alignments before calculating dN dS values
# usage: qsub generateKaKs_musclePal2nalCodeml.sh $subsetTag

# load necessary modules
module load bio/2.0

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

# set outputs directory name
outFolder=$inputsPath"/selectionTests"

# retrieve subset tag
subsetTag=$1

# set up final results file
#resultsFile=$outFolder"/kaksResults_"$subsetTag".csv"
#echo "geneID  t  S  N  dNdS  dN  dS" > "$resultsFile"

# set up gene lengths files
pepLengths=$outFolder"/geneLengths_"$subsetTag".pep.csv"
echo "geneID,reference,consensus" > "$pepLengths"
cdsLengths=$outFolder"/geneLengths_"$subsetTag".cds.csv"
echo "geneID,reference,consensus" > "$cdsLengths"

# status message
echo "Begining $subsetTag subset analysis..."

# make working directory
workingDir=$outFolder"/selectionTests_"$subsetTag
mkdir $workingDir

# move to directory of inputs for running codeml
cd $workingDir

# set fasta files of pep and cds seqs
tmpRefPep=$outFolder"/Pulex.pep.flt.fa."$subsetTag
tmpConPep=$outFolder"/Olympics.pep.flt.fa."$subsetTag
tmpRefNuc=$outFolder"/Pulex.cds.flt.fa."$subsetTag
tmpConNuc=$outFolder"/Olympics.cds.flt.fa."$subsetTag

# loop over all genes in the reference
while IFS= read -r line; do
	# retrieve gene tag
	gTag=$(echo "$line" | cut -f1 | sed 's/>//g')

	# Status message
	echo "Analyzing $gTag..."

	# set gene output paths
	gFile=$workingDir"/"$gTag"_Daphnia.pep.tmp.fa"
	gRefNuc=$workingDir"/"$gTag"_refNuc.cds.tmp.fa"
	gConNuc=$workingDir"/"$gTag"_conNuc.cds.tmp.fa"
	outAln=$workingDir"/"$gTag"_Daphnia_aligned.pep.tmp.fa"

	# calculate gene lengths for pep sequences
	pepRefSize=$(grep -w "^>$gTag" $tmpRefPep | cut -f2 | sed 's/\./X/g')
	pepConSize=$(grep -w "^>$gTag" $tmpConPep | cut -f2 | sed 's/\./X/g')
	echo $gTag","${#pepRefSize}","${#pepConSize} >> $pepLengths

	# calculate gene lengths for pep sequences
	cdsRefSize=$(grep -w "^>$gTag" $tmpRefNuc | cut -f2 | sed 's/\./X/g')
	cdsConSize=$(grep -w "^>$gTag" $tmpConNuc | cut -f2 | sed 's/\./X/g')
	echo $gTag","${#cdsRefSize}","${#cdsConSize} >> $cdsLengths

	# retrieve peptide sequences and convert back to multiline fasta format
	#echo ">Pulex_$gTag" > $gFile
	#grep -w "^>$gTag" $tmpRefPep | cut -f2 | sed 's/\./X/g' >> $gFile
	#echo ">Olympics_$gTag" >> $gFile
	#grep -w "^>$gTag" $tmpConPep | cut -f2 | sed 's/\./X/g' >> $gFile

	# retrieve nucleotide sequences and convert back to multiline fasta format
	#echo ">Pulex_$gTag" > $gRefNuc
	#grep -w "^>$gTag" $tmpRefNuc | cut -f2 >> $gRefNuc
	#echo ">Olympics_$gTag" > $gConNuc
	#grep -w "^>$gTag" $tmpConNuc | cut -f2 >> $gConNuc

	# prepare tree file
	#echo "(>Pulex_$gTag, >Olympics_$gTag);" > $workingDir"/"$gTag".tree"

	# prepare control file template from pal2nal
	# seqtype = 1 for codon alignments
	# runmode = -2 performs ML estimation of dS and dN in pairwise comparisons
	#cat $baseDir"/util/test.cnt" | sed "s/test\.codon/$gTag\.codon/g" | sed "s/test\.tree/$gTag\.tree/g" | sed "s/test\.codeml/$gTag\.codeml/g" > $workingDir"/"$gTag".cnt"

	# Status message
	#echo "Generating MSA for $gTag..."

	# create MSA using muscle
	# https://stackoverflow.com/questions/70769809/muscle-command-line-wrapper
	#muscle -in "$gFile" -out "$outAln"
	#muscle -align "$gFile" -output "$outAln"

	# status message
	#echo "Generating codon alignment for $gTag..."

	# usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
	#"$softwarePath"/pal2nal.pl "$outAln" "$gRefNuc" "$gConNuc" -output paml -nogap  >  $workingDir"/"$gTag".codon"
	#pal2nal.pl "$outAln" "$gRefNuc" "$gConNuc" -output paml -nogap  >  $workingDir"/"$gTag".codon"

	# status message
	#echo "Generating ka ks values for $gTag..."

	# run codeml to retrieve ka ks values
	# you can find the output of codeml in the .codeml file
	# Ks, Ka values are very end of the output file
	# enter yes in case of stop codon warning prompt
	# https://groups.google.com/g/pamlsoftware/c/HNx4O_YMHVA
	# https://stackoverflow.com/questions/11190004/automatically-press-enter-to-continue-in-bash
	#yes "" | codeml $gTag".cnt"

	# save ka ks values to final results file
	#kaks=$(tail -1 $gTag".codeml")
	#echo "$gTag  $kaks" >> $resultsFile

	# clean up
	#rm $workingDir/*
done < $tmpRefPep

# clean up
rm -r $workingDir
rm $tmpRefPep
rm $tmpConPep
rm $tmpRefNuc
rm $tmpConNuc

# status message
echo "$subsetTag subset analysis complete!"

