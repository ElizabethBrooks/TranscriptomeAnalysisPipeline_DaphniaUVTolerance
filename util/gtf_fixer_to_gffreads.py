#Script to replace end positions of a GTF file with the sequence size of a single line FASTA file.
#by Victor Gambarini
#usage: python gtf_fixer_to_gffread.py file.gtf file.fasta

import sys
file_gtf = open(sys.argv[1])
gtf = file_gtf.readlines()
file_fasta = open(sys.argv[2])
fasta = file_fasta.readlines()
output = open(str(sys.argv[1]).split("/")[-1].split(".")[0]+"_fixed.gtf", "w")
log = open(str(sys.argv[1]).split("/")[-1].split(".")[0]+"_log.txt", "w")
transcript = {}
for line in fasta:
	if line[0] == '>':
		header = line[1:].split(" ")[0].rstrip()
		seq_len = 0
	else:
		seq_len += len(str(line).rstrip())
		transcript[header] = seq_len 
for line in gtf:
	line_list = line.split("\t")
	header = line_list[0]
	transcript_end = int(line_list[4])
	transcript_size = int(transcript[header])
	if transcript_end >= transcript_size:
		line_list[4] = str(transcript_size-1)
		output.write("\t".join(line_list))
		log.write(header+"\t"+str(transcript_end-transcript_size)+"\n")
	else:
		output.write(line)
output.close()
log.close()