#!/bin/bash
#Script to run Rscripts that generate binned kMeans plots

#Plot merged data binned kMeans clustering
Rscript geneCounts_kMeans_binned.r $1 $2
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv/\.pdf/')
mv Rplots.pdf $outFile