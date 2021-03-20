#!/bin/bash

#Prepare multiline protein fasta to retrieve seqs
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpE05.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Y05_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpY05.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Y023_5_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpY023_5.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1R2_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpR2.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1PA_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpSierra.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/databases/PA42_v4.1_proteins/PA42.4.1.pep_longest.fasta | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA42_v4.1.fasta

#Retrieve selected protein sequences and convert back to multiline fasta format
grep "^>TRINITY_GG_5631_c0_g1_i1.p1" tmpE05.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5631_c0_g1_i1.p1.*/>E05_TRINITY_GG_5631_c0_g1_i1.p1/g' > gene4487_proteins_allDaphnia.fasta
grep "^>TRINITY_GG_5677_c0_g1_i1.p1" tmpY05.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5677_c0_g1_i1.p1.*/>Y05_TRINITY_GG_5677_c0_g1_i1.p1/g' >> gene4487_proteins_allDaphnia.fasta
grep "^>TRINITY_GG_5492_c0_g1_i1.p1" tmpY023_5.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5492_c0_g1_i1.p1.*/>Y023_5_TRINITY_GG_5492_c0_g1_i1.p1/g' >> gene4487_proteins_allDaphnia.fasta
grep "^>TRINITY_GG_5504_c1_g1_i1.p" tmpR2.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5504_c1_g1_i1.p.*/>R2_TRINITY_GG_5504_c1_g1_i1.p/g' >> gene4487_proteins_allDaphnia.fasta
grep "^>TRINITY_GG_3028_c0_g1_i1.p2" tmpPA.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_3028_c0_g1_i1.p2.*/>PA_TRINITY_GG_3028_c0_g1_i1.p2/g' >> gene4487_proteins_allDaphnia.fasta
grep "^>TRINITY_GG_5339_c1_g1_i1.p1" tmpSierra.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5339_c1_g1_i1.p1.*/>Sierra_TRINITY_GG_5339_c1_g1_i1.p1/g' >> gene4487_proteins_allDaphnia.fasta
grep "^>dp_gene4487-pep" tmpPA42_v4.1.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>dp_gene4487-pep.*/>PA42_v4.1_dp_gene4487-pep/g' >> gene4487_proteins_allDaphnia.fasta

#Clean up
rm tmp*.fasta

#Create MSA
muscle -in gene4487_proteins_allDaphnia.fasta -out gene4487_proteins_allDaphnia_aligned.fasta


#Prepare multiline transcript fasta to retrieve seqs
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98/cdhitEst | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpE05.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Y05_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98/cdhitEst | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpY05.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Y023_5_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98/cdhitEst | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpY023_5.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1R2_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98/cdhitEst | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpR2.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1PA_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98/cdhitEst | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98/cdhitEst | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpSierra.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/databases/PA42_v4.1_transcripts/PA42.4.1.cdna.longest.fasta | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA42_v4.1.fasta

#Retrieve selected nucelotide sequences and convert back to multiline fasta format
grep "^>TRINITY_GG_5631_c0_g1_i1" tmpE05.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5631_c0_g1_i1.*/>E05_TRINITY_GG_5631_c0_g1_i1/g' > gene4487_transcripts_allDaphnia.fasta
grep "^>TRINITY_GG_5677_c0_g1_i1" tmpY05.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5677_c0_g1_i1.*/>Y05_TRINITY_GG_5677_c0_g1_i1/g' >> gene4487_transcripts_allDaphnia.fasta
grep "^>TRINITY_GG_5492_c0_g1_i1" tmpY023_5.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5492_c0_g1_i1.*/>Y023_5_TRINITY_GG_5492_c0_g1_i1/g' >> gene4487_transcripts_allDaphnia.fasta
grep "^>TRINITY_GG_5504_c1_g1_i1" tmpR2.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5504_c1_g1_i1.*/>R2_TRINITY_GG_5504_c1_g1_i1/g' >> gene4487_transcripts_allDaphnia.fasta
grep "^>TRINITY_GG_3028_c0_g1_i1" tmpPA.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_3028_c0_g1_i1.*/>PA_TRINITY_GG_3028_c0_g1_i1/g' >> gene4487_transcripts_allDaphnia.fasta
grep "^>TRINITY_GG_5339_c1_g1_i1" tmpSierra.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5339_c1_g1_i1.*/>Sierra_TRINITY_GG_5339_c1_g1_i1/g' >> gene4487_transcripts_allDaphnia.fasta
grep "^>dp_gene4487-mRNA" tmpPA42_v4.1.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>dp_gene4487-mRNA.*/>PA42_v4.1_dp_gene4487-mRNA/g' >> gene4487_transcripts_allDaphnia.fasta

#Clean up
rm tmp*.fasta

#Create MSA
muscle -in gene4487_transcripts_allDaphnia.fasta -out gene4487_transcripts_allDaphnia_aligned.fasta


#Prepare multiline cds fasta to retrieve seqs
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpE05.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Y05_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpY05.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Y023_5_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpY023_5.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1R2_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpR2.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1PA_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/lodge3/ebrooks/rnaseq/sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyPA42_v4.1Trinity/decoded_transdecoder/Trinity.fasta.transdecoder.cds | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpSierra.fasta
cat /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/databases/PA42_v4.1_cds/PA42.4.1.cds_longest.fasta | sed ':a;N;$!ba;s/\n/NEWLINE/g' | sed 's/NEWLINE>/\n>/g' > tmpPA42_v4.1.fasta

#Retrieve selected coding sequences and convert back to multiline fasta format
grep "^>TRINITY_GG_5631_c0_g1_i1.p1" tmpE05.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5631_c0_g1_i1.p1.*/>E05_TRINITY_GG_5631_c0_g1_i1.p1/g' > gene4487_cds_allDaphnia.fasta
grep "^>TRINITY_GG_5677_c0_g1_i1.p1" tmpY05.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5677_c0_g1_i1.p1.*/>Y05_TRINITY_GG_5677_c0_g1_i1.p1/g' >> gene4487_cds_allDaphnia.fasta
grep "^>TRINITY_GG_5492_c0_g1_i1.p1" tmpY023_5.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5492_c0_g1_i1.p1.*/>Y023_5_TRINITY_GG_5492_c0_g1_i1.p1/g' >> gene4487_cds_allDaphnia.fasta
grep "^>TRINITY_GG_5504_c1_g1_i1.p" tmpR2.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5504_c1_g1_i1.p.*/>R2_TRINITY_GG_5504_c1_g1_i1.p/g' >> gene4487_cds_allDaphnia.fasta
grep "^>TRINITY_GG_3028_c0_g1_i1.p2" tmpPA.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_3028_c0_g1_i1.p2.*/>PA_TRINITY_GG_3028_c0_g1_i1.p2/g' >> gene4487_cds_allDaphnia.fasta
grep "^>TRINITY_GG_5339_c1_g1_i1.p1" tmpSierra.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>TRINITY_GG_5339_c1_g1_i1.p1.*/>Sierra_TRINITY_GG_5339_c1_g1_i1.p1/g' >> gene4487_cds_allDaphnia.fasta
grep "^>dp_gene4487-CDS" tmpPA42_v4.1.fasta | sed 's/NEWLINE/\n/g' | sed 's/^>dp_gene4487-CDS.*/>PA42_v4.1_dp_gene4487-CDS/g' >> gene4487_cds_allDaphnia.fasta

#Clean up
rm tmp*.fasta

#Create MSA
muscle -in gene4487_cds_allDaphnia.fasta -out gene4487_cds_allDaphnia_aligned.fasta