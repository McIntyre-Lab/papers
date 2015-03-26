#!/usr/bin/perl
#$ -cwd
# Marty McCrory, 2011-09-13
#
# formats a fastq file:
# trims the barcode (first 5 nt) off
# trims leading As off
# trims trailing Bs (in quality score) off
# trims off beginning and end primers
# eliminates any reads that are super-short
# prints what's left
# output is a FASTQ file

use strict;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('3:5:mbps:h', \%options);
if ($options{h}){Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis program \"fixes\" a fastq file.  The program can perform the following operations:\n";
	print $fh "\t-5 <int>: trim the first <int> characters off the read.\n";
	print $fh "\t-3 <int>: trim the last <int> characters off the read.\n";
	print $fh "\t-m: trims homopolymers off the front of the read (two or more consecutive nt of the same base call).\n";
	print $fh "\t-b: trims b-tails off the 3' end of the read (bases with quality score \"B\" using Phred 1.6 scoring).\n";
	print $fh "\t-p: trims primer/adapter off the 5' and 3' ends of the read.\n";
	print $fh "\t-s <int>: When all \"fixes\" are done, the script discards reads with length <int> or less.  Defaults to 10 if not specified.\n";
	print $fh "\nNote: all options are optional.  If a particular option is not specified, that particular \"fix\" will not be performed.\n";
}

my $barcode_3_prime_length = $options{"3"}; # set 3' barcode length to input
my $barcode_5_prime_length = $options{"5"}; # set 5' barcode length to input
my $short_reads_length = $options{s} ? $options{s} : 10; # set short read threshold to input, or default of 10

# $options{h} = trim homopolymers
# $options{b} = trim b-tails
# $options{p} = trim primer/adapter


my $num_mismatches_allowed = 1;
my $primer_end = "AGATCGGAAG"; # first 10 nt of adapter on end of sequence (AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG)
my $primer_begin = "CTTCCGATCT"; # last 10 nt of adapter on begin of sequence (ACACTCTTTCCCTACACGACGCTCTTCCGATCT)
my %read_length_counts;
my $num_trimmed_As;
my $num_trimmed_Bs;
my $num_trimmed_primer_end;
my $num_trimmed_primer_begin;
my $num_too_short_after_trimming;
my $num_all_reads;

my ($header, $readseq, $line3, $quality) = ("-1", "-1", "-1", "-1");

## added by mccrory, 2011-09-30
## adds flags for what trimmings took place
my %trim_flag;
#print STDERR "read_id,trimmed_homopolymers,trimmed_trailing_B_qual,trimmed_primer_3',trimmed_primer_5'\n";

while (<>) {
	chomp;
	if (/^@/) { # we're at the header line in the fastq file
		
		# trim the barcode off the front of the read and quality
		if ($options{"5"}) {
			for (my $i=0; $i<$barcode_5_prime_length; $i++) {
				$readseq =~ s/^.//;
				$quality =~ s/^.//;
			}
		}
			
		# trim 5' homopolymers, added by mccrory 2011-09-20
		if ($options{m}) {
			if ($readseq =~ /^AA+|^CC+|^GG+|^TT+/) {
				$readseq =~ s/$&//;
				$quality = substr ($quality, length $&);
				$trim_flag{$header}{homop} = 1;
			}
		}
			
		# trim the 3' barcode, added by mccrory 2011-09-13
		if ($options{"3"}) {
			for (my $i=0; $i<$barcode_3_prime_length; $i++) {
				$readseq =~ s/.$//;
				$quality =~ s/.$//;
			}
		}
		
		# if there are trailing Bs in the quality, trim the quality and read
		if ($options{b}) {
			if ($quality =~ /(B+)$/) {
				$quality =~ s/(B+)$//;
				$readseq = substr ($readseq, 0, (length $readseq) - (length $1));
				$num_trimmed_Bs++;
				$trim_flag{$header}{trail_B} = 1;1
			}
		}
		
		# remove END primer and everything after
		# remove BEGINNING primer and everything before
		if ($options{p}) {
			for (my $i=0; $i<(length($readseq) - (length $primer_end)); $i++) {
				my $temp_readseq = substr ($readseq, $i, length $primer_end);
				my $mm_end = mismatch_count($temp_readseq, $primer_end); # use hamming distance to calculate mismatches
				my $mm_begin = mismatch_count($temp_readseq, $primer_begin); # use hamming distance to calculate mismatches
				if ($mm_end <= $num_mismatches_allowed) {
					#print STDERR "match primer on END.  mm score: $mm_end.  temp_readseq: $temp_readseq.  readseq: $readseq, quality: $quality\n";
					$readseq =~ s/($temp_readseq.*)$//;
					$quality = substr ($quality, 0, (length $quality) - (length $1));
					$num_trimmed_primer_end++;
					$trim_flag{$header}{primer_end} = 1;
				}
			
				if ($mm_begin <= $num_mismatches_allowed) {
					#print STDERR "match primer on BEGINNING.  mm score: $mm_begin.  temp_readseq: $temp_readseq.  readseq: $readseq, quality: $quality\n";
					$readseq =~ s/^(.*?$temp_readseq)//;
					$quality = substr ($quality, length $1, (length $quality) - (length $1));
					$num_trimmed_primer_begin++;
					$trim_flag{$header}{primer_begin} = 1;
				}
			}
		}
			
		# print whatever's left
		if (length $readseq <= $short_reads_length) { $num_too_short_after_trimming++ }
		else { 
			print "$header\n$readseq\n$line3\n$quality\n" unless ($header eq "-1" || $readseq eq "" || $quality eq "-1" || $readseq eq "-1" || $line3 eq "-1" || length $readseq ne length $quality || length $quality > 200);
			## added by mccrory, 2011-09-30
			#print STDERR "$header,$trim_flag{$header}{homop},$trim_flag{$header}{trail_B},$trim_flag{$header}{primer_end},$trim_flag{$header}{primer_begin}\n"; 
		}
		
		$read_length_counts{length $readseq}++;
		$num_all_reads++;
		($header, $readseq, $line3, $quality) = ("-1", "-1", "-1", "-1");
		$header = $_;
	}
	elsif (/^[ACGNT]+$/) { # we're at the readseq line
		$readseq = $_;
	}
	elsif (/^\+/) { # we're at the third line
		$line3 = $_;
	}
	else { # we're at the phred quality line
		$quality = $_;
	}
}

#print STDERR "leading_As,$num_trimmed_As\n";
#print STDERR "trailing_Bs,$num_trimmed_Bs\n";
#print STDERR "primer_end,$num_trimmed_primer_end\n";
#print STDERR "primer_begin,$num_trimmed_primer_begin\n";
#print STDERR "reads_length_<=_4,$num_too_short_after_trimming\n";

my $sum_of_all_lengths = 0;
foreach my $len (sort {$a <=> $b} keys %read_length_counts) {
	$sum_of_all_lengths += ($len * $read_length_counts{$len});
}

my $avg_trimmed_read_length = $sum_of_all_lengths / $num_all_reads;

my $to_print = shift @ARGV;
$to_print =~ s/^.*?To_Florida\///;
print STDERR "fastq_file,all_reads,avg_trimmed_length,leading_As,trailing_Bs,primer_end,primer_begin,reads_length_lt5\n";


print STDERR "$to_print,$num_all_reads,$avg_trimmed_read_length,$num_trimmed_As,$num_trimmed_Bs,$num_trimmed_primer_end,$num_trimmed_primer_begin,$num_too_short_after_trimming\n";

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
