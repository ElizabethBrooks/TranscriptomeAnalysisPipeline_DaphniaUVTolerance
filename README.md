# TranscriptomeAnalysisPipeline_DaphniaUVTolerance
Bash shell scripts for analyzing short paired-end RNA sequence reads for Daphnia UV tolerance treatments.

To avoid uploading outputs to this repository, the outputs from these scripts be placed in separate folders in the same directory as the script folder. In other words, outputs from scripts in this repository folder (Directory/TranscriptomeAnalysisPipeline_DaphniaUVTolerance) will be placed in the directory that contains this repository folder (Directory/ScriptOutputFolder).

## RNA-seq Analysis Pipeline
![RNA-seq Analysis Pipeline](RNASeq_Workflow_DmelUV.png)

## Running Scripts on Servers
* To submit a trimming job to the queue: **qsub *SCRIPTNAME*.sh** 
* To submit an alignment or stats job to the queue: **qsub *SCRIPTNAME*.sh *FOLDERNAME_run0* ... *FOLDERNAME_runN***  
* To view the jobs you have submitted and corresponding task ID numbers: **qstat -u *USERNAME***
* To delete a job from the queue: **qdel *TASKIDNUMBER***

## Running Scripts Locally
* To compile the script before running: **chmod +x *SCRIPTNAME*.sh**
* To run a compiled trimming script: **./*SCRIPTNAME*.sh**
* To run a compiled alignment or stats script: **./*SCRIPTNAME*.sh *FOLDERNAME_run0* ... *FOLDERNAME_runN***

## Script Naming Format
Each script is named by the action and the primary software needed to perform the action.

## Running Any Script
Instructions for usage, with required inputs are given in the first few lines of each script. There is a text file with information about the inputs, and one for outputs in the **InputData** folder. This is where paths may be set for input and output directories.

## Pipeline Component Scripts
These are scripts that perform a single pipeline operation, and are located in the **Pipeline** directory.

### Quality Control
* QC_fastqc.sh
  * Output: **trimmed_run#/SAMPLENAME_fastqc_report.txt**  

### Adapter Trimming
* trimming_trimmomaticFastqc.sh
  * Output: **trimmed_run#**  

### Sequence Alignment
These scripts will accept any number of folders with reads trimmed using Trimmomatic. A minimum of one folder is expected as input.
* alignment_hisat2.sh
  * Input(s): ***trimmed_run0* ... *trimmed_runN***  
  * Output: **aligned_hisat2_run#**   

### Statistical Analysis
These scripts will accept a mix of folders with reads aligned using either HISAT2 or Tophat2. A minimum of one folder is expected as input.
* sorting_samtools.sh
  * Input(s): ***aligned_SOFTWARE_run0* ... *aligned_SOFTWARE_runN***  
  * Output: **sorted_samtools_run#** 
* counting_htseq.sh
  * Input(s): ***sorted_SOFTWARE_run0* ... *sorted_SOFTWARE_runN***  
  * Output: **counted_htseq_run#** 

## Pipeline Stage Scripts
These are scripts that perform all operations necessary for a stage of the pipeline.

### Adapter Trimming with Quality Control
* trimmingQC_trimmomaticFastqc.sh
  * Output: **trimmed_run#**  

## Optimized Scripts
These are scripts that have been optimized for running on ND CRC servers using distributed resources.

## Workflow Summary ##
1. Quality control check a sample with [FastQC][2] to identify the correct adapter library encoding of illumina pipeline used and corresponding phred.
   * If *Encoding = Illumina 1.5*, then the phred score is 64  
   * If *Encoding = Illumina 1.9*, then the phred score is 33  
2. Perform adapter trimming with [Trimmomatic][3] for paired-end data with two specified input files, and resulting in 4 output files. Output files consist of 2 files for the paired output where both reads survived the processing, and 2 for corresponding unpaired output where a read survived, but the partner read did not. Adapter trimming is achieved by the
   1. Removal of adapters: *ILLUMINACLIP:/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10*
   2. Removal of leading low quality bases with a score below 3: *LEADING:3*
   3. Removal of trailing low quality bases with a score below 3: *TRAILING:3*
   4. Scanning of reads with a 4-base wide sliding window and cutting when the average quality per base drops below 12: *SLIDINGWINDOW:4:15*
   5. Dropping of reads below 36 bases long: *MINLEN:36*
   6. Cutting of specified number of bases from the start of the read: *HEADCROP:13*
3. Quality control check the trimmed paired reads to determine if “reads are good enough” to proceed.
4. Map trimmed reads using a reference genome to perform sequence alignment with the [HISAT2][5] or [Tophat2][6] packages to
   1. Check the mapping efficiency of each job
   2. Prepare reads for sorting and counting
5. Assemble transcripts and quantify samples for differential expression analysis using [Cufflinks][7] or [HTSeq-count][8], depending on chosen statistical analysis package.
6. Perform statistical analysis by generating read counts with the Tuxedo or [EdgeR][9] pipelines to
   1. Statistically find differences in expression levels
   2. Generate an annotation based on the mapped reads
   3. Perform differential gene expression analysis on the mapped reads

## Required Software ##
* [FastQC][10]: A quality control tool for high throughput raw sequence data. It generates quality reports for NGS data and gives pass/fail results for the following checks: Per base sequence quality, Per sequence quality scores, Per base sequence content, Per base GC content, Per sequence GC content, Per base N content, Sequence length distribution, Sequence duplication levels, Overrepresented sequences, Kmer content. It also has a Graphic User Interface.
* [Trimmomatic][11]: A flexible read trimming tool for Illumina NGS data. It can trim adapter sequences, remove low-quality reads and bases.
* [HISAT2][12]: A fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome, and exome sequencing data) against the general human population (as well as against a single reference genome). The algorithm is based on HISAT and Bowtie2; uses a graph FM index (GFM) to index the genome before read mapping.
* [Tophat2][13]: A spliced read mapper for RNA-Seq. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons.
* [Bowtie2][14]: An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. Bowtie2 first extracts "seed" substrings in reads, aligns seeds in an ungapped way, and then performs extension in a gapped way.
* [Cufflinks][15]: It assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples. Assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples. It can be used in the pipeline with a protocol paper.
* [Cuffdiff][16]: Differential analysis of gene regulation at transcript resolution with RNA-seq. An algorithm that estimates expression at transcript-level resolution and controls for variability evident across replicate libraries.
* [Samtools][17]: Utilities for the Sequence Alignment/Map (SAM) format. SAMtools has multiple commands for processing SAM/BAM files. The sub-command "SAMtools-flagstat" can be used to print statistics for SAM/BAM files using the FLAG field.
* [HTSeq-count][18]: A package to count mapped reads for genomic features. It counts mapped reads for genomic features.
* [EdgeR][19]: Empirical Analysis of Digital Gene Expression Data. It performs differential expression analysis using read counts. It uses raw count data; implements a range of statistical methodology based on the negative binomial distributions, including empirical Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests.

  [1]: https://files.osf.io/v1/resources/twvc5/providers/osfstorage/5d000f49fea9230019808e67?mode=render
  [2]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
  [3]: http://www.usadellab.org/cms/?page=trimmomatic
  [4]: http://www.htslib.org/doc/#manual-pages
  [5]: https://ccb.jhu.edu/software/hisat2/manual.shtml#running-hisat2
  [6]: https://ccb.jhu.edu/software/tophat/index.shtml
  [7]: http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html
  [8]: https://htseq.readthedocs.io/en/release_0.11.1/counting.html
  [9]: https://bioconductor.org/packages/release/bioc/html/edgeR.html
  [10]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
  [11]: http://www.usadellab.org/cms/?page=trimmomatic
  [12]: https://ccb.jhu.edu/software/hisat2/manual.shtml#running-hisat2
  [13]: https://ccb.jhu.edu/software/tophat/index.shtml
  [14]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
  [15]: http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html
  [16]: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/
  [17]: http://www.htslib.org/doc/#manual-pages
  [18]: https://htseq.readthedocs.io/en/release_0.11.1/counting.html
  [19]: https://bioconductor.org/packages/release/bioc/html/edgeR.html
