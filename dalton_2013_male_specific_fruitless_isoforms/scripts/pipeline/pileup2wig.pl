#! /usr/bin/perl

#This script converts a SAMtools pileup file to a WIG file

use strict;
use warnings;

die "usage: pileup2WIG_jmf.pl <PILEUP> <TRACK_NAME> > <OUTPUT>\n" unless $#ARGV+1>1;

my $name = $ARGV[1];
print "track type=wiggle_0 name=\"$name\" visibility=full color=0,0,0\n";

open (PILEUP, "$ARGV[0]") or die "Could not open PILEUP: $!\n";

my $last_Chrom;
my $last_ID  = "null";

while (<PILEUP>) {
	my @cols = split(/\t/, $_);
	my $chrom = $cols[0];
	my $pos = $cols[1];
	my $depth = $cols[7];
	
	$chrom =~ s/chr//g;
	# remove the chrom string matching if not using human genome
	if ($last_ID eq $chrom) { # && $chrom =~ /(^\d+$)|(^X)|(^Y)/) {
        print "$pos\t$depth\n";
	}
	
	if ($last_ID ne $chrom) {
		$last_ID = $chrom;
		# remove the chrom string matching if not using human genome
		#if ($chrom =~ /(^\d+$)|(^X)|(^Y)/) {
		#	print "variableStep\tchrom=chr$chrom\n$pos\t$depth\n";
		#}
		#else {
			print "variableStep\tchrom=$chrom\n$pos\t$depth\n";
		#}
	}
	
}
	
close PILEUP;
