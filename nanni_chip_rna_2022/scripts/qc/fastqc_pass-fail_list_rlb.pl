#!/usr/bin/perl 

use strict;
use warnings ;


my @array;
my $BC=$ARGV[0];
my $file=$ARGV[1];
my $basic_stats;
my $base_qual;
my $seq_qual;
my $base_content;
my $base_GC;
my $seq_GC;
my $base_N;
my $seq_length;
my $seq_dup;
my $seq_overrep;
my $kmer;
my $count;

open (IN,"<","$file");
while (<IN>) {
	
	my @array=split(/\t/);
	if ($.==1) {
	$basic_stats=$array[0];
	}
	if ($.==2) {
	$base_qual=$array[0];
	}
	if ($.==3) {
        $seq_qual=$array[0];
	}
	if ($.==4) {
        $base_content=$array[0];
	}
	if ($.==5) {
        $base_GC=$array[0];
	}
	if ($.==6) {
	$seq_GC=$array[0];
	}
	if ($.==7) {
        $base_N=$array[0];
	}
	if ($.==8) {
        $seq_length=$array[0];
	}
	if ($.==9) {
        $seq_dup=$array[0];
	}
	if ($.==10) { 
        $seq_overrep=$array[0];
	}
	if ($.==11) {
        $kmer=$array[0];
	}

	if ($array[0] =~ /PASS/) {
	$count++;
	}
}

my $percent_pass=substr(($count/11)*100, 0, 5) ;

print join ",", $BC, $basic_stats, $base_qual, $seq_qual, $base_content, $base_GC, $seq_GC, $base_N, $seq_length, $seq_dup, $seq_overrep, $kmer, $percent_pass;

print "\n";

 
