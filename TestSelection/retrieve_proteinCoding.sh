#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveProteins_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub retrieve_proteinCoding.sh

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeFeatures://g")
#genomeFeatures="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff"

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" $baseDir"/InputData/outputPaths.txt" | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics"

# retrieve genome reference absolute path for alignment
refPath=$(grep "genomeReference" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")
#refPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# make outputs directory name
outFolder=$inputsPath"/features_gffread"

# retrieve input bam file type
type="filteredMapQ"

# retrieve file name of reference
refTag=$(basename $refPath)

# set paths for protein coding sequence lists
geneList=$outFolder"/"$refTag"_proteinCoding_genes.txt"
transList=$outFolder"/"$refTag"_proteinCoding_map.txt"

# set reference pep and cds paths
inRefPep=$outFolder"/"$refTag"_longest.pep.fa"
inRefNuc=$outFolder"/"$refTag"_longest.cds.fa"

# set consensus pep and cds paths
inConPep=$outFolder"/"$type"_consensus_longest.pep.fa"
inConNuc=$outFolder"/"$type"_consensus_longest.cds.fa"

# set reference multiline pep fasta to retrieve seqs
tmpRefPep=$outFolder"/Pulex.pep.tmp.fa"
fltRefPep=$outFolder"/Pulex.pep.flt.fa"

# set input consensus multiline pep fasta
tmpConPep=$outFolder"/Olympics.pep.tmp.fa"
fltConPep=$outFolder"/Olympics.pep.flt.fa"

# set reference multiline cds fasta to retrieve seqs
tmpRefNuc=$outFolder"/Pulex.cds.tmp.fa"
fltRefNuc=$outFolder"/Pulex.cds.flt.fa"

# set input consensus multiline cds fasta
tmpConNuc=$outFolder"/Olympics.cds.tmp.fa"
fltConNuc=$outFolder"/Olympics.cds.flt.fa"

# pre-clean up
rm $transList
rm $fltRefPep
rm $fltConPep
rm $fltRefNuc
rm $fltConNuc

# create singleline fasta of seqs
cat $inRefPep | sed 's/\ gene=.*/\t/g' | tr -d '\n' | sed 's/>/\n>/g' > $tmpRefPep
cat $inConPep | sed 's/\ gene=.*/\t/g' | tr -d '\n' | sed 's/>/\n>/g' > $tmpConPep
cat $inRefNuc | sed 's/\ gene=.*/\t/g' | tr -d '\n' | sed 's/>/\n>/g' > $tmpRefNuc
cat $inConNuc | sed 's/\ gene=.*/\t/g' | tr -d '\n' | sed 's/>/\n>/g' > $tmpConNuc

# create list of protein coding sequence gene names
cat $genomeFeatures | grep "biotype" | grep "protein_coding" | cut -d ";" -f 5 | sed 's/=/-/g' > $geneList

# loop over each gene name
while IFS= read -r line; do
	# status message
	echo "Processing $line ..."
	# create list of protein coding sequence transcript names, and grab the first listed transcript
	transName=$(cat $genomeFeatures | grep "$line" | grep "mRNA" | head -1 | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2)
	# add transcript name to the gene to transcript map file
	echo "$transName $line" >> $transList
	# reference pep fasta
	echo -en ">"$line"\t" >> $fltRefPep
	grep -w "$transName" $tmpRefPep | cut -f 2 >> $fltRefPep
	# consensus pep fasta
	echo -en ">"$line"\t" >> $fltConPep
	grep -w "$transName" $tmpConPep | cut -f 2 >> $fltConPep
	# reference cds fasta
	echo -en ">"$line"\t" >> $fltRefNuc
	grep -w "$transName" $tmpRefNuc | cut -f 2 >> $fltRefNuc
	# consensus cds fasta
	echo -en ">"$line"\t" >> $fltConNuc
	grep -w "$transName" $tmpConNuc | cut -f 2 >> $fltRefNuc
done < $geneList

# clean up
#rm $geneList
#rm $tmpRefPep
#rm $tmpConPep
#rm $tmpRefNuc
#rm $tmpConNuc

# status message
echo "Analysis complete!"
