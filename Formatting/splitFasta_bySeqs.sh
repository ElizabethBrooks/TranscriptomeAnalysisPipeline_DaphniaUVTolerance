#!/usr/bin/env Rscript

chunkNum=$2
seqNum=$(grep ">" "$1" | wc -l)
chunkSize=$((($seqNum($chunkNum-1))/$chunkNum)) #round up

