# TranscriptomeAnalysisPipeline_DaphniaUVTolerance
Bash shell scripts for analyzing short RNA sequence reads for Daphnia UV tolerance treatments.

These scripts are designed to be run from within a script folder, which should be placed in the directory with the .fq.gz files to be processed.

## RNA-seq Analysis Pipeline
![RNA-seq Analysis Pipeline](RNASeq_Workflow_DmelUV.png)

## Naming
Each script is named by the action and the necessary software.

## Pipeline Component Scripts
These are scripts that perform a single pipeline operation.

*Quality Control*
1. QC_fastqc.sh

*Adapter Trimming*
1. trimming_trimmomaticFastqc.sh

*Sequence Alignment*
1. alignment_hisat2.sh
2. alignment_tophat2.sh

## Pipeline Stage Scripts
These are scripts that perform all operations necessary for a stage of the pipeline.

*Adapter Trimming with Quality Control*
1. trimmingQC_trimmomaticFastqc.sh
