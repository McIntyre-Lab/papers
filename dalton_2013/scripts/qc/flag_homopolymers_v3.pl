#!/usr/bin/perl 
#===============================================================================
#
#         FILE: flag_homopolymers_v3.pl
#
#        USAGE: ./flag_homopolymers_v3.pl  <INPUT> > <OUTPUT>
#
#  DESCRIPTION: We are looking at MOSAIK and determining the best settings to
#  apply. One issue we are having is that we are finding reads that have no
#  mismatches, but have bad mosaik quality scores (0). This script is going to
#  flag the un-aligned reads from the alignment pipeline and flag them if they
#  have a homopolymer.
#
#       AUTHOR: Justin Fear (JMF), 
#  INSTITUTION: University of Florida
#      VERSION: 3.0
#      CREATED: 12/21/2011 08:39:01 AM
#     REVISION: 12/22/2011
#               - Edited to count homopolymers any where in the read too, so we
#                 can look at ambig reads.
#===============================================================================

use strict;
use warnings;

my $flag1 = 0; # Flag to show that I have hit a new header /^@.*/ line
my $flag2 = 0; # Flag to show that I have picked up seq information and am ready to print
my $homoflag = 0; # OUTPUT: Flag that I have a homopolymer
my $homoflagend = 0; # OUTPUT: Flag that I have a homopolymer
my $header; # OUTPUT: All of the header information except @
my $seq; # OUTPUT: All of the sequence information

print "read_id,sequence,flag_homopolymer,flag_homopolymer_end\n";
while(<>){
    chomp;
    if(m/^@(.*)$/){
        $flag1 = 1;
        $header = $1;
    }
    elsif ($flag1 == 1){
        $seq = $_;
        $flag1 = 0;
        $flag2 = 1;
        if (m/A{8,}$|C{8,}$|T{8,}$|G{8,}$/g){
            $homoflagend = 1;
        }
        if (m/A{8,}|C{8,}|T{8,}|G{8,}/g){
            $homoflag = 1;
        }
        
    }
    if ($flag2 == 1){
        print "$header,$seq,$homoflag,$homoflagend\n";
        $homoflag = 0;
        $homoflagend = 0;
        $flag2 = 0;
    }
}
