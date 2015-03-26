#! /usr/bin/perl
# Victor Amin 2010

use strict;
use warnings;

use FindBin '$Bin';
use lib "$Bin/modules";
use QueryFASTA;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('r:f:h', \%options);

if ($options{h} || !$options{r}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis program generates a FASTA file from BED.\nSTDIN/STDOUT\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-r [reference FASTA] [required]\n";
	print $fh "\t-f [flanking sequence] [default: 0]\n";
	print $fh "\n";
	exit;
}

my $flank = $options{f} ? $options{f} : 0;

my $qfa = new QueryFASTA($options{r});

print STDERR "\nConverting BED...\n";
my $count;
while (<>) {
	if (/^track\s/) {next}
	chomp;
	my ($chrom,$chrom_start,$chrom_end,$ident,$score,$strand,$t_start,$t_end,$rgb,$block_count,$block_sizes,$block_starts) = split;
	$chrom_start++; # 0-based to 1-based
	my @sizes = $block_sizes ? split(',', $block_sizes) : ();
	my @starts = $block_starts ? split(',', $block_starts) : ();
	
	my $seq = '';

	if ($flank > 0) {$seq = $qfa->getSequence($chrom, $chrom_start-$flank, $chrom_start-1); $seq .= "\n";}
	if (@sizes > 0) {
		for my $i (0..$#sizes) {
			my $start = $chrom_start + $starts[$i];
			my $end = $chrom_start + $starts[$i] + $sizes[$i] - 1;
			
			$ident .= "|$chrom:$start-$end";
			$seq .= $qfa->getSequence($chrom, $start, $end);
			$seq .= "\n";
		}
	} else {
		$ident .= "|$chrom:$chrom_start-$chrom_end";
		$seq .= $qfa->getSequence($chrom, $chrom_start, $chrom_end);
		$seq .= "\n";
	}
	if ($flank > 0) {$seq .= $qfa->getSequence($chrom, $chrom_end+1, $chrom_end+$flank); $seq .= "\n";}
	
	if ($seq eq '0') {print STDERR "WARNING: Sequence $ident|$chrom:$chrom_start-$chrom_end not in provided reference."} else {print ">$ident\n$seq"}
	$count++;
}

print STDERR "\nConverted $count rows.\n\n";

