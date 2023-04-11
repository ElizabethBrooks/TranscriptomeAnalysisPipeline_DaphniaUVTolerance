#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateKaKs_jobOutput

# script to generate MSAs for each gene in the reference set of peptide sequences
# usage: bash generateKaKs_musclePal2nalCodeml.sh
# usage ex: bash generateKaKs_musclePal2nalCodeml.sh

#Set software path
#softwarePath=$(grep "pal2nal:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/pal2nal://g")

#Retrieve sorted reads input absolute path
#inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/"
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/selectionTests"
mkdir $outFolder

# set inputs folder
inputsPath=$inputsPath"/features_gffread"

#Retrieve genome reference absolute path for alignment
#refPath=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
refPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# retrieve file name of reference
refTag=$(basename $refPath)

# set reference pep path
refPath=$inputsPath"/"$refTag"_longest_pep.fa"

# set consensus pep path
inputsPath=$inputsPath"/"$type"_consensus_longest_pep.fa"

# set results file path
resultsFile=$outFolder"/kaksResults.csv"

#Name output file of inputs
inputOutFile=$outFolder"/generateMSA_summary.txt"

#Add version to output file of inputs
muscle --version > $inputOutFile
#pal2nal.pl --version >> $inputOutFile

# status message
echo "Generating MSA for protein sequences..."

#Prepare reference multiline pep fasta to retrieve seqs
tmpRef=$outFolder"/tmpPulex_pep.fa"
cat $refPath | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpRef

#Prepare input multiline pep fasta
tmpSample=$outFolder"/tmpOLYM_pep.fa"
cat $inputsPath | sed 's/$/NEWLINE/g' | tr -d '\n' | sed 's/NEWLINE>/\n>/g' > $tmpSample

#Loop over all genes in the reference
while IFS= read -r line; do
	#Retrieve selected peptide sequences and convert back to multiline fasta format
	gTag=$(echo $line | sed 's/NEWLINE/\n/g' | grep ">" | sed 's/>//g')
	gFile=$outFolder"/tmp_"$gTag"_Daphnia_pep.fa"
	grep "^>$gTag" $tmpRef | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>Pulex_$gTag/g" > "$gFile"
	
	# Status message
	echo "Generating MSA for $gTag..."

	#Prepare multiline pep fasta to retrieve seqs
	grep "^>$gTag" $tmpSample | sed 's/NEWLINE/\n/g' | sed "s/^>$gTag.*/>OLYM_$gTag/g" >> "$gFile"

	# set file paths
	outAln=$outFolder"/tmp_"$gTag"_Daphnia_aligned_pep.fa"
	inAln=$outFolder"/"$gTag"_Daphnia_pep.fa"

	# create MSA using muscle from biopython
	# https://stackoverflow.com/questions/70769809/muscle-command-line-wrapper
	muscle -align "$gFile" -output "$outAln"

	# replace X wildcard with stop codon *
	sed 's/X/\*/g'  "$outAln" > "$inAln"

	# clean up
	rm $gFile
	rm $outAln

	# prepare tree file
	echo "(>PA42_v4.1_$gTag, >Olympics_$gTag);" > "$outFolder"/"$gTag".tree

	#Prepare control file
	# seqtype = 1 for codon alignments
	# runmode = -2 performs ML estimation of dS and dN in pairwise comparisons
	cat ../util/test.cnt | sed "s/test\.codon/$gTag\.codon/g" | sed "s/test\.tree/$gTag\.tree/g" | sed "s/test\.codeml/$gTag\.codeml/g" > $outFolder"/"$gTag".cnt"

	# status message
	echo "Generating codon alignment for $gTag..."

	#Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
	#"$softwarePath"/pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  $outFolder"/"$gTag".codon"
	pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  $outFolder"/"$gTag".codon"

	#Move to directory of inputs for codeml
	#cd $outFolder

	# status message
	echo "Generating ka ks values for $gTag..."

	#Run codeml to retrieve ka ks values
	#You can find the output of codeml in the .codeml file
	#Ks, Ka values are very end of the output file
	codeml $outFolder"/"$gTag".cnt"

	#Save ka ks values to final results file
	kaks=$(tail -1 $gTag".codeml")
	echo "$gTag  $kaks" >> $resultsFile

	#Clean up
	[ -f "rst" ] && rm "rst"
	[ -f "rst1" ] && rm "rst1"
	[ -f "rub" ] && rm "rub"
	rm $gTag".codon"
	rm $gTag".tree"
	rm $gTag".cnt"
done < $tmpRef

# clean up
rm $tmpSample
rm $tmpRef

# status message
echo "MSA generated for protein sequences!"

