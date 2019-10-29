#!/bin/bash
#Script to run Rscripts that generate kmeans PCA plots

#Plot merged data kmeans PCA
Rscript geneCounts_kMeans_binned.r ../../GeneCounts_Merged/final_merged_counts_annotated_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_merged_kMeans.pdf

#Plot fullset data kmeans PCA
Rscript geneCounts_kMeans.r ../../GeneCounts_Merged/merged_counts_fullset_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_fullset_kMeans.pdf

#Plot subset data kmeans PCA
Rscript geneCounts_kMeans.r ../../GeneCounts_Merged/merged_counts_subset_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_subset_kMeans.pdf

#Plot legacy data kmeans PCA
Rscript geneCounts_kMeans.r ../../GeneCounts_Merged/merged_counts_legacy_transposed.csv
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_legacy_kMeans.pdf