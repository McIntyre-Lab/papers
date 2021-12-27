#!/usr/bin/perl 
#===============================================================================
#
#         FILE: macs_bed2sbed.pl
#
#        USAGE: ./macs_bed2sbed.pl  
#
#  DESCRIPTION: This file takes the bed and summit bed file from MACS and makes
#  a single bed file.
#
#      OPTIONS: 
#       AUTHOR: Justin Fear (JMF),jfear@ufl.edu 
# ORGANIZATION: University of Florida
#      VERSION: 1.0
#      CREATED: 12/19/2012 01:55:02 PM
#===============================================================================

use strict;
use warnings;
use Getopt::Long;

my $bed = '';
my $summit = '';

GetOptions('bed=s' => \$bed, 'summit=s'=> \$summit );

my %storage = ();

open(BED,"<",$bed);

while(<BED>){
    my($chr,$start,$end,$name,$rank) = split;
    $storage{$name}{'chr'} = $chr;
    $storage{$name}{'start'} = $start;
    $storage{$name}{'end'} = $end;
    $storage{$name}{'summit'} = '';
}
close BED;

open(SUMMIT, "<", $summit);

while(<SUMMIT>){
    my($chr,$start,$end,$name,$rank) = split;
    $storage{$name}{'summit'} = $end;
}
close SUMMIT;

foreach my $key (keys %storage){
    print "$storage{$key}{'chr'}\t$storage{$key}{'start'}\t$storage{$key}{'end'}\t$storage{$key}{'summit'}\t$key\n";
}
