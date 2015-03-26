#! /usr/bin/perl
# Victor Amin 2010

use strict;
use warnings;

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('b:p:h', \%options);

if ($options{h} || !$options{b}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis script splits alignments to junctions into two separate alignments to chromosomes. STDIN/STDOUT, counts to STDERR.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-b <file> junctions.bed [required]\n";
	print $fh "\t-p <int> junctions padding [default: read length (autodetected)]\n";
	print $fh "\n\n";
	exit;
}

my $padding = $options{p} ? $options{p} : -1;
open JUNCBED, $options{b} or die "Could not open junctions .bed file: $!\n";

my %juncs;
while (<JUNCBED>) {
	if (/^track/) {next}
	chomp;
	my @cols = split;
	
	my $chrom = $cols[0];
	my $chrom_start = $cols[1] + 1; # 0 -> 1 based
	my $junc_id = $cols[3];
	my @lengths = split(',', $cols[10]);
	my @starts = split(',', $cols[11]); # 0 based
	
	$juncs{$junc_id}{chrom} = $chrom;
	$juncs{$junc_id}{chrom_start} = $chrom_start;
	$juncs{$junc_id}{break_pos} = $lengths[0] + 1;
	$juncs{$junc_id}{second_exon_start} = $chrom_start + $starts[1];
	$juncs{$junc_id}{chrom_end} = $chrom_start + $starts[1] + $lengths[1] - 1;
}

my $count = 0;
my $nj = 0;
while (<>) {
	if (/^@/) {next}
	chomp;

	my @cols = split;
	
	my $query_name = $cols[0];
	my $flag = $cols[1];
	my $subject_name = $cols[2];
	my $pos = $cols[3];
	my $mapq = $cols[4];
	my $mrnm = $cols[6];
	my $mpos = $cols[7];
	my $isize = $cols[8];
	my $seq = $cols[9];
	my $qual = $cols[10];
	my @tags = @cols[11..$#cols];
	my $strand = '+';

	if ($flag == 4) {next}
	
	if ($padding == -1) {
		$cols[5] =~ /(\d+)M/;
		$padding = $1;
	}
	
	my $binflag = dec2bin($flag);
	my @flags;
	my $i = 0;
	while (length $binflag > 0) {
		$flags[2**$i] = chop $binflag;
		$i++;
	}
	
	if (exists $flags[16] && $flags[16] == 1) {$strand = "-"}
	
	$count++;

	my @s = split(/\|/, $subject_name);
	my $junc_id = $s[0];

	my $length = length($seq);
	my $chrom = $juncs{$junc_id}{chrom};
	my $chrom_pos = $juncs{$junc_id}{chrom_start} + $pos - $padding - 1;
	my $length_1 = $juncs{$junc_id}{break_pos} - $pos + $padding;
	my $length_2 = $length - $length_1;
	
	my $xsa = 0;
	for my $i (0..$#tags) {
		if ($tags[$i] =~ /XS:A:./) {
			$xsa = 1;
		}
	}
	if ($xsa == 0) {push @tags, "XS:A:$strand"}
	my $tags = join("\t", @tags);

	if ($length_1 <= 0) {
		$nj++;
		# after the break
		$chrom_pos = $juncs{$junc_id}{second_exon_start} + $pos - $juncs{$junc_id}{break_pos} - $padding;
		print "$query_name\t$flag\t$chrom\t$chrom_pos\t$mapq\t${length}M\t$mrnm\t$mpos\t$isize\t$seq\t$qual\t$tags\n";
		next;
	} elsif ($length_1 >= $length) {
		$nj++;
		# before the break
		print "$query_name\t$flag\t$chrom\t$chrom_pos\t$mapq\t${length}M\t$mrnm\t$mpos\t$isize\t$seq\t$qual\t$tags\n";
		next;
	}
	
	my $N = $juncs{$junc_id}{second_exon_start} - $chrom_pos - $length_1;
	print "$query_name\t$flag\t$chrom\t$chrom_pos\t$mapq\t${length_1}M${N}N${length_2}M\t$mrnm\t$mpos\t$isize\t$seq\t$qual\t$tags\n";
}

print STDERR "Aligned to junctions: $count\nAligned near junctions: $nj\n";


## FUNCITONS

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

