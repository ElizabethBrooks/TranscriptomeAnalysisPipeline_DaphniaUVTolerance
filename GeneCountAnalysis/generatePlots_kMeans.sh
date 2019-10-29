#!/bin/bash
#Script to run Rscripts that generate kMeans plots

#Plot merged data kMeans clustering
Rscript geneCounts_kMeans.r $1
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv/\.pdf/')
mv Rplots.pdf $outFile