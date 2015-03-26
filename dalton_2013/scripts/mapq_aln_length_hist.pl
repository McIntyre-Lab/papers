#! /usr/bin/perl

use strict;

my %aln_qualities;
print "aln_quality,aln_length,num_reads\n";
while (<>) {
	chomp;
	my ($read_id, $aln_quality, $aln_length) = split ","; 
	
	$aln_qualities{$aln_quality}{$aln_length}++;
}

for my $aln_quality (sort {$a <=> $b} keys %aln_qualities) {
	for my $aln_length (sort {$a <=> $b} keys %{$aln_qualities{$aln_quality}}) {
		print "$aln_quality,$aln_length,$aln_qualities{$aln_quality}{$aln_length}\n";
	}
}
