#!/bin/bash
#Script to convert a Genbank fna file to fasta format

cat /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fna | sed "s/^>.*scaffold/>scaffold/g" | sed "s/,.*//g" > /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta