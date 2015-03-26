#!/usr/bin/perl 
#===============================================================================
#
#        USAGE: ./identify_motifs_in_genes_v6.pl  
#
#  DESCRIPTION: This script takes a list of motif locations throughout the
#  geneome and compares it to the locations of genes of interest.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#        NOTES: ---
#       AUTHOR: Justin Fear (JMF), jfear@ufl.edu
# ORGANIZATION: University of Florida
#      VERSION: 1.0
#      CREATED: 08/21/2012 06:21:44 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use Getopt::Std;
#use Data::Dumper;

# Input Options
$Getopt::Std::STANDARD_HELP_VERSION = 1;
my %options;
getopts('m:g:u:d:h', \%options);
if ($options{h}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}

my $motif_file = $options{m};
my $gene_file = $options{g};

# number of bases upstream
my $ulength = $options{u};

# number of bases downstream
my $dlength = $options{d};

my %gene_hash = ();
my $count;
my $cpos;
my $ceval;
open(GENE,"<","$gene_file");
while(<GENE>){
    if($. == 1){next}
    chomp;
    my ($symbol,$fbgn,$chrom,$start,$end,$strand) = split(",",$_);
    $gene_hash{$chrom}{$fbgn}= {
                                  symbol    => $symbol,
                                  start     => $start,
                                  end       => $end,
                                  strand    => $strand,
                                  mstart    => $start - $ulength,
                                  mend      => $end + $dlength,
                                  count     => 0,
                                  position  => '',
                                  eval      => ''
                               };
}
close(GENE);

open(MOTIF,"<","$motif_file");
while(<MOTIF>){
    chomp;
    if($. == 1){next}
    my($chrom,$pos,$eval,$match) = split(",",$_);

    # Keep only perfect matches, if there is a space in the match string, skip
    if($match =~ / /){next}

    # Don't need to look at Uextra
    if($chrom =~ /Uextra/){next}

    if(!exists $gene_hash{$chrom}){print "ERROR: Chromosome Name ($chrom) Does not exist"; exit}

    foreach my $fbgn (keys $gene_hash{$chrom}){

        if($pos >= $gene_hash{$chrom}{$fbgn}{mstart} && $pos <= $gene_hash{$chrom}{$fbgn}{mend}){

            if($gene_hash{$chrom}{$fbgn}{count} != 0){

                $cpos  = $gene_hash{$chrom}{$fbgn}{position} . ";$pos";
                $ceval = $gene_hash{$chrom}{$fbgn}{eval} . ";$eval";
                $count = $gene_hash{$chrom}{$fbgn}{count} + 1;

            } else{
                $count = 1;
                $cpos  = $pos;
                $ceval  = $eval;
            }

            $gene_hash{$chrom}{$fbgn}{count}    = $count;
            $gene_hash{$chrom}{$fbgn}{position} = $cpos;
            $gene_hash{$chrom}{$fbgn}{eval}     = $ceval;
        }
    }

}
close(MOTIF);

print "primary_fbgn,symbol,chrom,start,end,strand,upstream_start,upstream_end,motif_count,motif_positions,motif_eval\n";
foreach my $chrom (keys %gene_hash){
    foreach my $fbgn (keys $gene_hash{$chrom}){
        my $symbol = $gene_hash{$chrom}{$fbgn}{symbol};
        my $start  = $gene_hash{$chrom}{$fbgn}{start};
        my $mstart = $gene_hash{$chrom}{$fbgn}{mstart};
        my $end    = $gene_hash{$chrom}{$fbgn}{end};
        my $mend   = $gene_hash{$chrom}{$fbgn}{mend};
        my $pos    = $gene_hash{$chrom}{$fbgn}{position};
        my $eval   = $gene_hash{$chrom}{$fbgn}{eval};
        my $count  = $gene_hash{$chrom}{$fbgn}{count};
        my $strand = $gene_hash{$chrom}{$fbgn}{strand};

        print "$fbgn,$symbol,$chrom,$start,$end,$strand,$mstart,$mend,$count,$pos,$eval\n";
    }

}

#-------------------------------------------------------------------------------
#  Subroutines
#-------------------------------------------------------------------------------

sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "This script takes a list of motif locations throughout the geneome and compares it to the locations of genes of interest.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-m Motif file from parsed MAST xml output.\n";
	print $fh "\t-g gene list information\n";
	print $fh "\n";
	exit;
}

