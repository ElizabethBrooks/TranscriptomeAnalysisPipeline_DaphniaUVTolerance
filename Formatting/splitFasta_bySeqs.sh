#!/bin/bash
#Script to split an input fasta file input appropriate number of sequence
# chunks based on the input max chunk size (in MB)
#Usage: splitFasta_bySeqs.sh fastaFile maxMB
#Usage Ex: splitFasta_bySeqs.sh ../../Daphnia_pulex_PA42_v3.0.fasta 50

#Retrieve input fasta absolute path
inputsPath=$(dirname $1)
#Determine the number of pieces to split the file into
chunkSize=$2 #Max number of MB for each chunk
fastaSize=$(grep ">" "$1" | wc -c)
fastaSize=$(($fastaSize*0.000001)) #Convert bytes to MB
chunkNum=$((($fastaSize+($chunkSize-1))/$chunkSize)) #Round up
#Determine the number of sequences to include in each chunk
seqsNum=$(grep ">" "$1" | wc -1)
seqsNum=$((($seqsNum+($chunkNum-1))/$chunkNum)) #Round up
#Prepare a temporary fasta file for splitting with cut
# by replacing all newlines with a marker
sed ':a;N;$!ba;s/\n/NEWLINE/g' > tmp.fasta
#Retrieve the input fatsa file name
inFileNoEx=$(basename $1)
inFileNoEx=$(echo $inFileNoEx | sed 's/\.fasta//')
#Loop through the input fasta file and split by sequence chunks
fileNum=0
startPos=0
endPos=$seqsNum
for (( i=0; i<=$chunkNum; i++ )); do
	echo "Creating sequnce chunk $i file..."
	#Output current sequence chunk to file
	fileNum=$(($i+1))
	cut -d "<" -f $startPos-$endPos tmp.fasta > "$inputsPath"/"$inFileNoEx"_"$fileNum".fasta
	echo "Sequnce chunk $i file created!"
	#Calculate the sequence chunk window shift
	window=$(($seqsNum*$i))
	startPos=$(($startPos+$window))
	endPos=$(($endPos+$window))
done
#Clean up
rm tmp.fasta