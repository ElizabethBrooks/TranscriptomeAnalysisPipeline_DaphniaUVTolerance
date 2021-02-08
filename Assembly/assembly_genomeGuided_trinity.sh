#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assemblyGenome_trinity_jobOutput
#$ -pe smp 8
#Script to perform genome-guided Trinity de novo transcriptome assembly
#Usage: qsub assembly_genomeGuided_trinity.sh sortedFolder genotype maxIntronLength genome
#Usage Ex: qsub assembly_genomeGuided_trinity.sh sortedCoordinate_samtoolsHisat2_run2 E05 14239 PA42_v3.0
#Usage Ex: qsub assembly_genomeGuided_trinity.sh sortedCoordinate_samtoolsHisat2_run2 E05 23554 PA42_v4.1

#Required modules for ND CRC servers
module load bio
#module load bio/samtools
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No arguments supplied... exiting"
   	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != sortedCoordinate* ]]; then
	echo "ERROR: The "$1" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$1" == *"Hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$1" == *"Tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The sorted "$1" folder of bam files were not found... exiting"
	exit 1
fi
#Retrieve aligned reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#Retrieve assembly outputs absolute path
outputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
#Create output directory
outputFolder="$outputsPath"/"$1""$2"_assembly"$4"Trinity
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to outputs directory
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1""$2"_assembly"$4"Trinity_summary.txt
#Merge and re-coordinate sort the set of bam files
readFiles=$(echo "$inputsPath"/"$1"/*_"$2"_*/*.bam)
echo "Beginning merging..."
samtools merge -@ 8 merged.bam $readFiles
echo "Merging complete! Beginning sorting..."
samtools sort -@ 8 -o sorted.bam merged.bam
echo "Sorting complete!"
rm merged.bam
#Run Trinity on coordinate-sorted bam files using 8 threads, and a maximum intron
# length that makes most sense given your targeted organism
#Reduce output with --full_cleanup
echo "Beginning assembly of $1 reads for $2 data..."
Trinity --genome_guided_bam sorted.bam --genome_guided_max_intron "$3" --max_memory 50G --CPU 8 --output "$outputFolder"
echo "Assembly complete!"
rm sorted.bam
#Add run inputs to output summary file
echo "samtools merge --threads 8 merged.bam "$readFiles > "$inputOutFile"
echo "samtools sort -@ 8 -o sorted.bam merged.bam" >> "$inputOutFile"
echo "Trinity --genome_guided_bam sorted.bam --genome_guided_max_intron" "$3" "--max_memory 50G --CPU 8 -output" "$outputFolder" >> "$inputOutFile"
#Remove selected extra files generated by Trinity
rm -r "$outputFolder"/Dir_sorted.bam.minC1.gff
rm "$outputFolder"/Dir_sorted*
rm "$outputFolder"/sorted.bam*
rm "$outputFolder"/partitions.ok
rm "$outputFolder"/read_files.list*
rm "$outputFolder"/trinity_GG.cmds*
#Re-name selected output files
mv "$outputFolder"/Trinity-GG.fasta.gene_trans_map "$outputFolder"/Trinity.fasta.gene_trans_map
mv "$outputFolder"/Trinity-GG.fasta "$outputFolder"/Trinity.fasta
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"