#! /usr/bin/perl
# Marty McCrory 2011-09-19

# input: .bed file, .pileup file



use strict;
no strict "refs";
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('b:p:o:', \%options);

open BED, $options{b};
open PILEUP, $options{p};

my $OUTDIR = $options{o};

my %regions;
while (<BED>) {
	chomp;
	my ($chrom, $start, $end, $id, $score, $strand) = split;
	for (my $i=$start; $i<$end; $i++) {
		$regions{$chrom}{$i} = $id;
		#print STDERR "$chrom,$i,$regions{$chrom}{$i}\n";
	}
	open ("$id", ">", "$options{o}/$id.consensus.pileup") or die "could not open file $options{o}/$id.pileup: $!";
}
close BED;

while (<PILEUP>) {
	my @cols = split;
	if ($regions{$cols[0]}{$cols[1]} ne "") { # we're at a pileup location that corresponds to one of the regions in the .bed file
		my $id = $regions{$cols[0]}{$cols[1]};
		#print STDERR "$id,$cols[0],$cols[1],$regions{$cols[0]}{$cols[1]}\n";
		print {$id} "$_";
	}
}
