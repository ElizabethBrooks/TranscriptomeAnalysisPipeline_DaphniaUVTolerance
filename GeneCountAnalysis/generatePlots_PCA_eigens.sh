#!/bin/bash
#Script to run Rscripts that generate PCA plots

#Plot merged data PCA with eigenvectors
Rscript geneCounts_PCA_eigens_binned.r ../../GeneCounts_Merged/final_merged_counts_annotated_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_merged_eigens_PCA.pdf

#Plot fullset data PCA with eigenvectors
Rscript geneCounts_PCA_eigens.r ../../GeneCounts_Merged/merged_counts_fullset_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_fullset_eigens_PCA.pdf

#Plot subset data PCA with eigenvectors
Rscript geneCounts_PCA_eigens.r ../../GeneCounts_Merged/merged_counts_subset_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_subset_eigens_PCA.pdf

#Plot legacy data PCA with eigenvectors
Rscript geneCounts_PCA_eigens.r ../../GeneCounts_Merged/merged_counts_legacy_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_legacy_eigens_PCA.pdf