#!/bin/bash
#Script to run Rscripts that generate PCA plots

#Plot merged data PCA
Rscript geneCounts_PCA_merged.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_merged_PCA.pdf
#Plot merged data PCA with eigenvectors
Rscript geneCounts_PCA_merged_eigens.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_merged_eigens_PCA.pdf

#Plot subset data PCA
Rscript geneCounts_PCA_fullset.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_fullset_PCA.pdf
#Plot merged data PCA with eigenvectors
Rscript geneCounts_PCA_fullset_eigens.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_fullset_eigens_PCA.pdf

#Plot merged data PCA
Rscript geneCounts_PCA_subset.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_subset_PCA.pdf
#Plot merged data PCA with eigenvectors
Rscript geneCounts_PCA_subset_eigens.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_subset_eigens_PCA.pdf

#Plot subset data PCA
Rscript geneCounts_PCA_legacy.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_legacy_PCA.pdf
#Plot merged data PCA with eigenvectors
Rscript geneCounts_PCA_legacy_eigens.r
#Rename produced plot
mv Rplots.pdf ../../GeneCounts_Merged/geneCounts_legacy_eigens_PCA.pdf