## Workflow Command Summary ##

*ND CRC access and data management*

 1. ssh userName@crcfe01.crc.nd.edu
 2. cd /afs/crc.nd.edu/group/hoth
 3. mkdir directoryName
 4. cp -avr /afs/crc.nd.edu/group/pfrenderlab/devries/melanica/rnaseq/dmelUV directoryName
 5. cd directoryName/dmelUV
 
*Quality Control*
 1. module load bio
 2. fastqc sample.fq.gz --extract
 3. if grep -iF "Encoding	Illumina 1.5" sample_fastqc/fastqc_data.txt; then score=64
 4. if grep -iF "Encoding	Illumina 1.9" sample_fastqc/fastqc_data.txt; then score=33

*Adapter trimming*
 1. module load bio/trimmomatic/0.32
 2. mkdir trimmed
 3. trimmomatic PE -"phred"$score sample_1 sample_2 trimmed/sample_pForward.fq.gz trimmed/sample_uForward.fq.gz trimmed/sample_pReverse.fq.gz trimmed/sample_uReverse.fq.gz ILLUMINACLIP:/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13

*Quality Control*
 1. fastqc sample_pForward.fq.gz --extract
 2. grep -iF "pattern"
 3. trimmed/sample_pForward_fastqc/fastqc_data.txt > trimmed/pattern.txt
 4. grep -iF "pattern" trimmed/sample_pReverse_fastqc/fastqc_data.txt > trimmed/pattern.txt
 
*Sequence alignment*
- Tophat2 Package
     1. mkdir aligned
     2. bowtie2-build /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned
     3. tophat2 -G genomeFilePath -o aligned/out/sample indexedGenomeFilePath trimmed/sample_pForward.fq.gz trimmed/sample_pReverse.fq.gz
     4. grep -iF "pattern" aligned/out/sample/align_summary.txt > aligned/sample_summaryReport.txt

- HISAT2 Package
     1. module load bio/hisat2/2.1.0
     2. mkdir aligned
     3. hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa Daphnia_pulex.allmasked
    hisat2 -q -x indexedGenomeFilePath -1 trimmed/sample_pForward.fq.gz -2 trimmed/sample_pReverse.fq.gz -S aligned/out/sample --aligned/sample_summaryReport.txt

*Statistical analysis*
- Tuxedo Pipeline
     1. module load bio/cufflinks/2.2.1
     2. cuffdiff -o stats genomeFilePath aligned/sample1/accepted_hits.bam, aligned/sample2/accepted_hits.bam, …, aligned/sampleN/accepted_hits.bam

- EdgeR Pipeline
     1. module load bio/python/2.7.14
     2. module load bio/htseq/0.11.2
     3. samtools sort -T /tmp/sample/accepted_hits_sorted.bam -o aligned/sample/accepted_hits_sorted.bam aligned/sample/accepted_hits.bam
     4. htseq-count -s no -m union -t gene -i trID aligned/sample/accepted_hits_sorted.bam -i genomeFilePath > sample.counts

## Workflow Command Reference ##

|Command|Description|
|------|----------|
|[fastqc][6] filePath|Check quality of reads with FastQC|
|[fastqc][6] filePath --extract|Extract files after running FastQC|
|[trimmomatic][8] PE -phred64 inputForward inputReverse pairedOutputForward.fq.gz unpairedOutputForward.fq.gz pairedOutputReverse>.fq.gz unpairedOutputReverse.fq.gz step1 … stepN|Trim paired end reads with Trimmomatic and specified phred|
|[bowtie2][9]-build genomeFilePath outputFolderPath|Align sequences to reference genome|
|[hisat2][10] [options]* -x hisat2-idx {-1 forwardReadPath -2 reverseReadPath} [-S alignedReadsOutputPath]|Map trimmed paired reads with HISAT2|
|[tophat2][11] -G genomeFilePath -o outputDirectoryName forwardReadPath reverseReadPath|Map trimmed reads with Tophat2 and a reference genome|
|[cufflinks][1] [options] alignedReadsPath|Generate annotations based on trimmed reads with Cufflinks|
|[cuffdiff][2] [options] transcripts.gtf sample1 sample2 … sampleN|Perform differential expression analysis with Cuffdiff|
|[samtools][4] sort [options] BAMFilePath|Sort BAM files with Samtools|
|[htseq][5]-count -s no -m method -t featureType -i id bamFilePath gffFilePath > output.counts|Generate counts for BAM file with HTSeq|

## General Command Reference ##

|Command|Description|
|------|----------|
|ssh userName@crcfe01.crc.nd.edu|Secure connection to ND CRC servers|
|ssh -X userName@crcfe01.crc.nd.edu|Enable graphics to open html files|
|ls|List files in current directory|
|cd pathTo|Change directories|
|cd ..|Move up one directory|
|mkdir folderName|Make a new directory|
|cp -avr pathFrom pathTo|Recursively copy files in directory|
|scp pathFrom pathTo|Copy files and directories to another server|
|rm filePath|Remove files|
|rm -rf folderPath|Remove directory and contents without prompts|
|grep -iF “[pattern]” inputFileName > outputFileName|Find the line in a file with the specified pattern and output the line to a new file|
|gzip -cd filePath1 filePath2 … filePathN|Concatenate compressed files|
|gzip > mergedFilePath|Recompress merged files|
|gunzip filePath|Decompress a compressed file|
|cat filePath1 filePath2 … filePathN > mergedFilePath|Merge list of files and write to a new file|
|module load bio|Load bio module with software for workflow|
|module list|List available modules|
|qsub filePath|Submit job|
|qstat -u userName|Display job status|
|qdel jobNumber|Delete specified job|
|firefox filePath|Open html file in browser|
|exit|Close server connection|

  [1]: http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html
  [2]: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html
  [3]: https://bioconductor.org/packages/release/bioc/html/edgeR.html
  [4]: http://www.htslib.org/doc/#manual-pages
  [5]: https://htseq.readthedocs.io/en/release_0.11.1/counting.html
  [6]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
  [8]: http://www.usadellab.org/cms/?page=trimmomatic
  [9]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
  [10]: https://ccb.jhu.edu/software/hisat2/manual.shtml#running-hisat2
  [11]: https://ccb.jhu.edu/software/tophat/index.shtml
