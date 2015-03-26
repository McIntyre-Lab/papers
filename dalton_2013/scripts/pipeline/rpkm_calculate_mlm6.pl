#! /usr/bin/perl
# Victor Amin 2010
# Modified Marty McCrory, April 2011

use strict;
#use warnings;

# added from rpkm_calc_mlm5.pl:
# outputs sample name also, does not output header

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('b:g:p:s:i:h', \%options);

if ($options{h} || (!$options{b} && !$options{g}) || !$options{p} || !$options{s}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis script calculates RPKM for any input. STDIN/STDOUT, counts to STDERR.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-b <file> input BED\n";
	print $fh "\t-g <file> input GTF\n";
	print $fh "\t-p <file> consensus pileup [required]\n";
	print $fh "\t-s <file> alignments SAM [required]\n";
	print $fh "\t-i <name> sample_id [required]\n";
	print $fh "\nPileup is from SAMtools, (consensus) variety. Output is a CSV with columns ID,MAPPED_READS,READS_IN_EXON,COVERAGE_IN_EXON,EXON_LENGTH,APN,RPKM\n";
	exit;
}

open SAM, "<$options{s}" or die "Could not open sam: $!\n";
open PILEUP, "<$options{p}" or die "Could not open pileup:$!\n";
if ($options{b}) {open INFILE, "<$options{b}" or die "Could not open ranges: $!\n"}
elsif ($options{g}) {open INFILE, "<$options{g}" or die "Could not open GTF: $!\n"}

my $sample_id = $options{i} ? $options{i} : "sample";
#print "ID,SAMPLE_ID,MAPPED_READS,READS_IN_EXON,COVERAGE_IN_EXON,EXON_LENGTH,APN,RPKM\n";
print STDERR "Reading SAM file...\n";
my $mapped_reads = 0;
my $read_length = 0;
while (<SAM>) {
	if (/^@/) {next}
	my ($read_id, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq) = split;
	if ($flag != 4) {$mapped_reads++}
	if ($. < 10000) { # autodetect read length
		my $length = length $seq;
		if ($read_length < $length) {$read_length = $length}
	}
}
close SAM;

my %pileup_depth;

my $id = '';
while (<PILEUP>) { # uses consensus pileup WITH LAST TWO COLUMNS TRIMMED
	my ($chrom, $pos, $ref, $con, $qual, $snp, $qual2, $depth, $junk1, $junk2) = split;
	$pileup_depth{$chrom}{$pos} += $depth;
}
close PILEUP;

print STDERR "Reading BED/GTF file and calculating coverage...\n";
while (<INFILE>) {
		chomp;
		my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes, $id);
		if ($options{g}) { ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = split; }
		elsif ($options{b}) { ($seqname, $start, $end, $id) = split; }
		
		my ($exon_length, $coverage, $reads_in_exon, $rpkm, $apn);

		$exon_length = ($end - $start + 1);
		for (my $i=$start; $i<=$end; $i++) { $coverage += $pileup_depth{$seqname}{$i};	}
		$reads_in_exon = $coverage/$read_length; 
		$rpkm = (1000000000 * $reads_in_exon) / ($mapped_reads * $exon_length); # from Mortazavi et al., Nature Methods 2008, supplement
		$apn = $coverage / $exon_length;
		print "$id,$sample_id,$mapped_reads,$reads_in_exon,$coverage,$exon_length,$apn,$rpkm\n";
}
close INFILE;

