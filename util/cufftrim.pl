#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 cufftrim.pl genome.fa.fai transcripts.gtf > transcripts.trim.gtf
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $fai=shift(@ARGV) || die ("$usage\nSamtools fasta index required!\n");
my $gtf=shift(@ARGV) || die ("$usage\nCufflinks' transcripts file required!\n");
my %slen;

open(FAI, $fai) || die("Error opening $fai\n");
open(GTF, $gtf) || die("Error opening $gtf\n");
while (<FAI>) {
 chomp;
 next if m/^#/;
 my @t=split(/\t/);
 die("Parsing error: $t[1] length for sequence $t[0] ?\n") 
    unless $t[1]>10;
 $slen{$t[0]}=int($t[1]);
}
close(FAI);

# now process the GTF file
my ($tcount, $ecount); #trim counts
while (<GTF>) {
 my @t=split(/\t/);
 my $l=$slen{$t[0]} || die("Error: contig $t[0] not found in the fasta index!\n");
 if ($t[4]>$l) {
   $t[4]=$l;
   $t[2] eq 'exon' ? ++$ecount : ++$tcount;
 }
 print join("\t", @t);
}
close(GTF);
print STDERR "$tcount transcripts, $ecount exons adjusted.\n";
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }