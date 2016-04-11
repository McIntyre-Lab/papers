#!/usr/bin/perl
# Marty McCrory, 2011-11-07

# input: weight matrix of sequences (columns in A T C G order)
# output: simulated sequences that follow this weight matrix, in fasta format
# Revision:
#           - Marty had the column order wrong it should be (A T G C) corrected JMF
#           - Marty actually had the column order correct, but he had the
#           column names wrong. It is actually (A T C G). Corrected JMF

use strict;

my %sequences;
while (<>) {
	chomp;
	next if /^A/;
	my ($a_percent, $t_percent, $c_percent, $g_percent) = split;
	my $total = 0;
	for my $i (0..$a_percent-1) { $sequences{$i} .= "A" }
	$total += $a_percent;
	for my $i ($total..$total+$t_percent-1) { $sequences{$i} .= "T" }
	$total += $t_percent;
	for my $i ($total..$total+$c_percent-1) { $sequences{$i} .= "C" }
	$total += $c_percent;
	for my $i ($total..$total+$g_percent-1) { $sequences{$i} .= "G" }
	$total += $g_percent;
}

# Meme fails to work if the sequence is less than 8 characters long. This loop
# adds an extra base on.
while (length $sequences{99} < 8) { # if the sequence isn't long enough
	for my $i (0..24) { $sequences{$i} .= "A" }
	for my $i (25..49) { $sequences{$i} .= "T" }
	for my $i (50..74) { $sequences{$i} .= "C" }
	for my $i (75..99) { $sequences{$i} .= "G" }
}

for my $seq (sort keys %sequences) {
	print ">$seq\n$sequences{$seq}\n";
}
