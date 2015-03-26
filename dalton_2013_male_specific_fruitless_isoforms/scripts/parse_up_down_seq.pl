#!/usr/bin/perl 
#===============================================================================
#
#        USAGE: ./parse_up_down_seq.pl  
#
#  DESCRIPTION: The ultimate upstream or downstream parsing script.
#
#      OPTIONS: 
#               -u extract upstream sequence. 1 extracts upstream sequence, 0 does not. [default: 1]
#               -d extract downstream sequence. 1 extracts downstream sequence, 0 does not. [default: 0]
#               -i Input FASTA file
#               -c Coordinate file
#               -p Plus Strand Indicator [default: +]
#               -m Minus Strand Indicator [default: -]
#               -l Length of region to grab [default: 1000]
#
# REQUIREMENTS: Getopt::Std
#               Bio::SeqIO
#
#        NOTES: ---
#       AUTHOR: Justin Fear (JMF), jfear@ufl.edu
# ORGANIZATION: University of Florida
#      VERSION: 1.0
#      CREATED: 08/20/2012 03:27:20 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;

# Input Options
$Getopt::Std::STANDARD_HELP_VERSION = 1;
my %options;
getopts('u:d:i:c:p:m:l:h', \%options);
if ($options{h}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}











#-------------------------------------------------------------------------------
#   Subroutines
#-------------------------------------------------------------------------------

sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis is the ultimate upstream/downstream parsing script. It takes in a FASTA file and gene coordiantes and pulls our the upstream and or downstream sequence into a new FASTA file.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-u extract upstream sequence. 1 extracts upstream sequence, 0 does not. [default: 1]\n";
	print $fh "\t-d extract downstream sequence. 1 extracts downstream sequence, 0 does not. [default: 0]\n";
	print $fh "\t-i Input FASTA file\n";
	print $fh "\t-c Coordinate file\n";
	print $fh "\t-p Plus Strand Indicator [default: +]\n";
	print $fh "\t-m Minus Strand Indicator [default: -]\n";
	print $fh "\t-l Length of region to grab [default: 1000]\n";
	print $fh "\n";
	exit;
}

sub import_coords{
}


sub plus_extract {
    if(@_ != 5){die "Malformed id2location";}
    my($ts,$chrom,$start,$end) = @_;
    my $ustart = $start - ($length + 1);
    my $upseq = substr($seq_hash{$chrom},$ustart,$length);

    print OUT ">$ts|$chrom $start $end plus_strand|$length bases upstream\n$upseq\n\n";
}

sub minus_extract {
    if(@_ != 5){die "Malformed id2location";}
    my($ts,$chrom,$start,$end) = @_;
    my $ustart = $end;
    my $upseq = substr($seq_hash{$chrom},$ustart,$length);
    my $revup = reverse($upseq);
    $revup =~ tr/ATCGatcg/TAGCtagc/ ;

    #open(OUT,">>","/home/jfear/sorghum/reports/1kb_upstream/upstream_seqs_v2.fa");
    print OUT ">$ts|$chrom $start $end minus_strand|$length bases upstream\n$revup\n\n";
    #close(OUT);
}
