## FastQC Notes ##
- **Base Sequence Quality:** The majority of the read should be green, and trail off into the yellow quality score range at the end.
- **Per Tile Sequence Quality:** This allows you to identify machine issues since it allows you to look at the quality scores from each tile across all of your bases to see if there was a loss in quality.
- **Per Sequence Quality Scores:** Only a small percentage of the total sequences should have universally poor quality.
- **Per Base Sequence Content:** The lines in this plot should run parallel with each other, and the relative amount of each base should reflect the overall amount of these bases in your genome. Typically, a warning or error will be given for reads produced with libraries that use random hexamers and those which were fragmented using transposases.
- **Per Sequence GC Content:** The distributions should overlap for the majority of mean GC content. An unusually shaped distribution could indicate a contaminated library or some other kinds of biased subset.
- **Per Base N Content:** If the proportion of Ns in the sequence rises above a few percent, it suggests that the analysis pipeline was unable to interpret the data well enough to make valid base calls. A very low proportion of Ns may be expected near the end of sequences.
Sequence Length Distribution: For some sequencing platforms it is entirely normal to have different read lengths so warnings can be ignored.
- **Sequence Duplication Levels:** A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication is more likely to indicate some kind of enrichment bias.
- **Overrepresented Sequences:** Finding that a single sequence is very overrepresented in the set either means that it is highly biologically significant, or indicates that the library is contaminated, or not as diverse as you expected.
- **Adapter Content:** This can find a number of different sources of bias in the library which can include the presence of read-through adapter sequences building up on the end of your sequences. However, finding the presence of any overrepresented sequences in your library will cause the Kmer plot to be dominated by the Kmers these sequences contain.
- **Kmer Content:** The Kmer module starts from the assumption that any small fragment of sequence should not have a positional bias in its appearance within a diverse library, and any Kmers with positionally biased enrichment are reported.

----------

## Manual ##
- [About][2]
- [Opening a Sequence File][3]
- [Evaluating Results][4]
- [Saving a Report][5]
- [Basic Statistics][6]
- [Per Base Sequence Quality][7]
- [Per Sequence Quality Scores][8]
- [Per Base Sequence Content][9]
- [Per Sequence GC Content][10]
- [Per Base N Content][11]
- [Sequence Length Distribution][12]
- [Duplicate Sequences][13]
- [Overrepresented Sequences][14]
- [Adapter Content][15]
- [Kmer Content][16]
- [Per Tile Sequence Quality][17]

----------

## Example Reports ##
- [Good Illumina data][18]
- [Bad Illumina data][19]
- [Adapter dimer contaminated run][20]
- [Small RNA with read-through adapter][21]
- [Reduced Representation BS-Seq][22]
- [PacBio][23]
- [454][24]


  [1]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  [2]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/1%20Introduction/1.1%20What%20is%20FastQC.html
  [3]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/2%20Basic%20Operations/2.1%20Opening%20a%20sequence%20file.html
  [4]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/2%20Basic%20Operations/2.2%20Evaluating%20Results.html
  [5]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/2%20Basic%20Operations/2.3%20Saving%20a%20Report.html
  [6]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/1%20Basic%20Statistics.html
  [7]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html
  [8]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html
  [9]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html
  [10]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html
  [11]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html
  [12]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html
  [13]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html
  [14]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html
  [15]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html
  [16]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html
  [17]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html
  [18]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html
  [19]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html
  [20]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html
  [21]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/small_rna_fastqc.html
  [22]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/RRBS_fastqc.html
  [23]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/pacbio_srr075104_fastqc.html
  [24]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/454_SRR073599_fastqc.html
