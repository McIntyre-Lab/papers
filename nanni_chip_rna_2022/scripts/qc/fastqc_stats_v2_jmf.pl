#!/usr/bin/perl

use strict;
use warnings ;

my $header='';
my $NAME=$ARGV[0];
my $OUTDIR=$ARGV[1];
my $O1 = "$OUTDIR"."/"."$NAME"."_overrep.txt";
my $O2 = "$OUTDIR"."/"."$NAME"."_kmer.txt";


open (OUT1, ">>","$O1");
open (OUT2, ">>","$O2");
open (IN,"<","fastqc_data.txt");

while (<IN>) {

    if (/^>>END/) {
        next;
    }

    if (/^>>(.*)\t(.*)/) {
	$header=$1;
    }

    if ($header eq "Overrepresented sequences") {
        print OUT1 $_;
    }

    if ($header eq "Kmer Content"){
        print OUT2 $_;
    }

}

