#!/usr/bin/perl 
#===============================================================================
#
#        USAGE: ./parse_last_sam_v2.pl  <SAM FILE> <AMBIG> <UNIQ>
# #  DESCRIPTION: This script takes a SAM file from LAST output and parses reads
#  taking the alignment with the best alignment score (AS:i:##). If there are two
#  alignments with the same score then it outputs these as ambiguous.
#
#       AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  INSTITUTION: University of Florida
#      VERSION: 1.0
#      CREATED: 02/24/2012 02:31:13 PM
#     REVISION: Changed the logic so that I keep uniq aln, and ambig aln with
#     best score. 
#===============================================================================

use strict;
use warnings;

my $input=$ARGV[0];
my $ambig=$ARGV[1];
my $uniq=$ARGV[2];

open(INPUT,"<","$input");
open(AMBIG,">","$ambig");
open(UNIQ,">","$uniq");

my %holder;
my $ambig_count = 0;
my $uniq_count = 0;

while(<INPUT>){
    chomp;
    if (/^(\S*)\t.*AS:i:(\d+)/){

        my $header = $1;
        my $score = $2;

        push(@{$holder{$header}},{'score' => $score, 'line' => $_})

    }
}

foreach my $key (keys %holder){

    my @sorted = sort {$b->{'score'} <=> $a->{'score'}} @{$holder{$key}};
    my $len = @sorted;

    if ($len > 1){
        if($sorted[0]{'score'} == $sorted[1]{'score'}){

            $ambig_count++;

            for (my $i = 0; $i < $len; $i++){
                if($sorted[0]{'score'} > $sorted[$i]{'score'}){
                    last;
                } else {
                    print AMBIG "$sorted[$i]{'line'}\n";
                }
            }

        } else {
            print UNIQ "$sorted[0]{'line'}\n";
            $uniq_count++;
        }
    } else {
        print UNIQ "$sorted[0]{'line'}\n";
        $uniq_count++;
    }
}

print "LAST Ambiguously Aligned Reads: $ambig_count\n";
print "LAST Uniquely Aligned Reads: $uniq_count\n";
