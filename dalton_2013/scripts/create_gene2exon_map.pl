#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  create_gene2exon_map.pl
#
#        USAGE:  ./create_gene2exon_map.pl  
#
#  DESCRIPTION:  A script to clean up output from
#  /share/mclab/useful_dmel_data/flybase526/fb526_si_fusion_2_gene_id.sas7bdat
#  for use in the R_wiggle pipeline
#
#===============================================================================

use strict;
use warnings;

my $flybase526 = "/home/Justin/Documents/R_wiggle/data/Fb526_si_fusion_2_gene_id.txt";
my $outfile = "/home/Justin/Documents/R_wiggles/data/gene2exon_flybase526.txt";

open(FH, "$flybase526") or die;
open(OUT,">","$outfile") or die;

print OUT "gene_symbol\tfusion_id\ttranscript_id\tchormosome\tstart_exon\tend_exon\n";
while(<FH>){
    chomp;
    my @line = split("\t");
    print OUT "$line[15]\t$line[0]\t$line[11]\t$line[2]\t$line[8]\t$line[9]\n";

}
close FH;
close OUT;

