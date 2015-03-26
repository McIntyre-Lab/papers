#!/usr/bin/perl
# Marty McCrory, 2011-11-10

# input: mast output xml file
# output: .csv file containing all hits

my $current_chrom;

print "chrom,motif_hit_position,motif_hit_eval,match\n";
while (<>) {
	chomp;
	if (/^\s+\<sequence id=".*name="(.*?)"/) { # <sequence id=   tag
		$current_chrom = $1
	} 
	if (/^\s+\<hit pos="(.*?)".*pvalue="(.*?)".*match="(.*?)"/) {
		my $hit_pos = $1;
		my $p_value = $2;
		my $match = $3;
		print join ",", $current_chrom, $hit_pos, $p_value, $match;
		print "\n";
	}
}
