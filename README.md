# TranscriptomeAnalysisPipeline_DaphniaUVTolerance
Bash shell scripts for analyzing short RNA sequence reads for Daphnia UV tolerance treatments.

These scripts are designed to be run from within a script folder, which should be placed in the directory with the .fq.gz files to be processed.

## RNA-seq Analysis Pipeline
![RNA-seq Analysis Pipeline](RNASeq_Workflow_DmelUV.png)

## Running Scripts
- To submit a job to the queue:
$ qsub SCRIPTNAME.sh
- To view the jobs you have submitted and corresponding task ID numbers:
$ qstat -u USERNAME
- To delete a job from the queue:
$ qdel TASKIDNUMBER

## Naming
Each script is named by the action and the necessary software.

## Pipeline Component Scripts
These are scripts that perform a single pipeline operation.

*Quality Control*
- QC_fastqc.sh

*Adapter Trimming*
- trimming_trimmomaticFastqc.sh

*Sequence Alignment*
- alignment_hisat2.sh
- alignment_tophat2.sh

## Pipeline Stage Scripts
These are scripts that perform all operations necessary for a stage of the pipeline.

*Adapter Trimming with Quality Control*
- trimmingQC_trimmomaticFastqc.sh

## Optimized Scripts
These are scripts that have been optimized for running on ND CRC servers using distributed resources.
