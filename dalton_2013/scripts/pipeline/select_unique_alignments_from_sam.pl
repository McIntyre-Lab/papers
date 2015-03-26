#!/usr/bin/perl
# Martin McCrory, December 2011

# input: sam file
# output: sam file containing only reads that DID align uniquely

use strict;

my %reads;

while (<>) {
	my @cols = split;
	$reads{$cols[0]}{count}++;
	$reads{$cols[0]}{line} .= $_;
	$reads{$cols[0]}{aln} = $cols[1];
}

for my $read_id (keys %reads) {
	if ($reads{$read_id}{count}==1 && $reads{$read_id}{aln} ne "4") { print "$reads{$read_id}{line}"; }
}
