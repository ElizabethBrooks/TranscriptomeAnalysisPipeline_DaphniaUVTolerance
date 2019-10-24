#!/bin/bash
#Script to run Rscripts that generate kmeans PCA plots

#Plot merged data kmeans PCA
Rscript geneCounts_kMeans_merged.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_merged_kMeans.pdf

#Plot fullset data kmeans PCA
Rscript geneCounts_kMeans_fullset.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_fullset_kMeans.pdf

#Plot subset data kmeans PCA
Rscript geneCounts_kMeans_subset.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_subset_kMeans.pdf

#Plot legacy data kmeans PCA
Rscript geneCounts_kMeans_legacy.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_legacy_kMeans.pdf