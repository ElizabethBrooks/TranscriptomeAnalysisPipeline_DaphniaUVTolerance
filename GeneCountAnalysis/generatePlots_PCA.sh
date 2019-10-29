#!/bin/bash
#Script to run Rscripts that generate PCA plots

#Plot merged data PCA
Rscript geneCounts_PCA_binned.r ../../GeneCounts_Merged/final_merged_counts_annotated_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_merged_PCA.pdf

#Plot fullset data PCA
Rscript geneCounts_PCA.r ../../GeneCounts_Merged/merged_counts_fullset_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_fullset_PCA.pdf

#Plot subset data PCA
Rscript geneCounts_PCA.r ../../GeneCounts_Merged/merged_counts_subset_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_subset_PCA.pdf

#Plot legacy data PCA
Rscript geneCounts_PCA.r ../../GeneCounts_Merged/merged_counts_legacy_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_legacy_PCA.pdf