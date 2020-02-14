#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_trinity_jobOutput
#$ -pe smp 8
#Script to generate a trinity assembly for input transcripts
#Usage: qsub assembly_trinity.sh
#Usage Ex: qsub assembly_trinity.sh

#fastqc reads
#trim
#use on trimmed files!!!

#move to trimmed directory

#Load required modules for ND CRC servers
module load bio/2.0

#TO DO: update for flexible inputs and add directed outputs (assembled_trinity_run1)
Trinity --seqType fq --max_memory 50G --left *_pForward.fq --right *_pReverse.fq --CPU 8

#do fastqc again!