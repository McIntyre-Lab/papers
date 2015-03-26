#! /usr/bin/perl
# Marty McCrory 2011-09-19

# input: all alignment log files (bowtie first, then mosaik)

# output: date, lane, total_reads, aln_nth_pass, percent_aln_nth_pass, ambig_nth_pass, percent_ambig_nth_pass;

use strict;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('d:l:', \%options);

my $date = $options{d};
my $lane = $options{l};
my $total_reads;
my $num_aligned;
my $num_reads_into_last;
print "$date,$lane";
while (<>) {
	# bowtie parsing
	if (/reads processed: (\d+)/ && $. == 1) { 
		$total_reads = $1;
		print ",$1" 
	}
	if (/reads with at least one reported alignment: (\d+)/ && $. % 5 == 2) { 
		my $percent = $1/$total_reads * 100;
		$num_aligned += $1;
		print ",$1,$percent" 
	}
	if (/reads with alignments suppressed due to -m: (\d+)/ && $. % 5 == 4) { 
		my $percent = $1/$total_reads * 100;
		$num_aligned += $1;
		print ",$1,$percent" 
	}
	if (/# reads processed: (\d+)/ && ($. == 6 || $. == 11)) {
		print ",$1"
	}
	if (/# reads that failed to align: (\d+)/ && $. == 13) {
		$num_reads_into_last = $1;
	}
	# last parsing
	if (/^Number of reads uniquely aligned by LAST: (\d+)/) { 
		my $percent = $1/$total_reads * 100;
		$num_aligned += $1;
		print ",$num_reads_into_last,$1,$percent" 
	}
	
	# all original reads parsing
	if (/^,(\d+),/ && $. == 18) {
		$total_reads = $1;
		print ",$1"
	}
}
my $num_aligned_percent = $num_aligned/$total_reads * 100;
print ",$num_aligned_percent\n";

