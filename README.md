# TranscriptomeAnalysisPipeline_DaphniaUVTolerance
Bash shell scripts for analyzing short RNA sequence reads for Daphnia UV tolerance treatments.

# Naming
Each script is named by the action and the necessary software.

# Pipeline Component Scripts
These are scripts that perform a single pipeline operation.

*Quality Control*
1. QC_fastqc.sh

*Adapter Trimming*
1. trimming_trimmomaticFastqc.sh

*Sequence Alignment*
1. alignment_hisat2.sh
2. alignment_tophat2.sh

# Pipeline Stage Scripts
These are scripts that perform all operations necessary for an entire stage of the pipeline.

*Adapter Trimming with Quality Control*
1. trimmingQC_trimmomaticFastqc.sh
