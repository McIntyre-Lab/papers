#!/usr/bin/perl 
#===============================================================================
#
#         FILE: parse_gene_pileup.pl
#
#        USAGE: ./parse_gene_pileup.pl  
#
#  DESCRIPTION: This script will ...
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 10/03/2012 05:58:03 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('c:f:g:h', \%options);

if ($options{h} || !$options{c}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}

my $COORD = $options{c};
my $FILE = $options{f};
my $GENE = $options{g};

my($ID, $SEC_ID, $CHROM, $START, $END, $STRAND);

open(COR,"<","$COORD");
while(<COR>){
    chomp;
    my ($tID,) = split(',',$_); 
    if($tID eq $GENE){
        ($ID, $SEC_ID, $CHROM, $START, $END, $STRAND) = split(',',$_); 
    }
}
close(COR);

open(IN,"<",$FILE);
while(<IN>){
    chomp;
    my($cCHROM,$POS,$COUNT) = split(/,/,$_);
    if( ($CHROM eq $cCHROM and $START <= $POS and $END >= $POS) or $.==1){
        print "$_\n";
    }
}
close(IN);


#-------------------------------------------------------------------------------
#  Subroutines
#-------------------------------------------------------------------------------

sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis program ...\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-c ID to coordinate file [required]\n";
	print $fh "\t-f File in the format CHROM,POS,COUNT [required]\n";
	print $fh "\t-g gene in the exact format as in the ID to coordinate file [required]\n";
	print $fh "\t-h Print these options\n";
	print $fh "\n";
	exit;
}

