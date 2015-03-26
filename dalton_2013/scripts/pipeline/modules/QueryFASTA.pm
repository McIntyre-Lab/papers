#! /usr/bin/perl
# Victor Amin 2010

package QueryFASTA;

use strict;
use warnings;

# This module allows querying of FASTA databases by coordinate or sequence
#
# Initializing:
#	my $qfa = new QueryFASTA($fasta) # $fasta is either the path to a FASTA file or a directory containing FASTA files
#
# Public Methods:
#	getHeaders: no params, returns a reference to an array of all FASTA headers
#	getLength: (headerName), returns length of sequence with that header
#	getNucleotide: (headerName, 1-basedCoord), returns the nucleotide at a position
#	getSequence: (headerName, 1-basedStartCoord, endCoord), returns sequence by coords
#	getSequenceSubstr: (headerName, 0-basedStartCoord, length), returns sequence, same syntax as Perl substr command
#	findSequence: (sequence), returns reference to hash of all exact matches for a sequence,
#											hash has format: (headerName => 1-basedStartCoord = endCoord)
#	freeSequence: (headerName), deletes a sequence from memory
#
#	All methods return 0 on failure
#

## CONSTRUCTOR/DESTRUCTOR ##
sub new
{
	my ($class, $fasta_location, $param_ref) = @_;
	my $self = {
		_mode => 'frontload',
		_fasta_location => $fasta_location,
		_messages => '1'
	};	# _mode: may one day implement on-disk access, but only if memory becomes an issue
		# _messages: 0 -- no logging to STDERR, 1 -- simple 'Loading' message, 2 -- list all FASTA records as they are loaded
	bless $self, $class;
	if ($param_ref) {
		$self->setParams($param_ref);
	}
	$self->_load();
	return $self;
}

## GETTERS/SETTERS ##
sub getParam
{
	my $self = shift;
	my ($param) = @_;
	return $self->{$param};
}

sub setParam
{
	my $self = shift;
	my ($param, $setting) = @_;	
	if ($param) {$self->{$param} = $setting} else {return 0}
	return $self->{$param};
}

sub setParams
{
	my $self = shift;
	my ($params_ref) = @_;
	if ($params_ref) {
		while(my($param,$setting) = each(%$params_ref)) {
			$self->{$param} = $setting;
		}
		return \(keys %$params_ref);
	}
	return 0;
}

## QUERY FASTA ##
sub getHeaders
{
	my $self = shift;
	my @headers = keys %{ $self->{_sequences} };
	return \@headers;
}

# $N = chrom name
# $S = start
# $E = end
# $L = length

sub getLength
{
	my $self = shift;
	my ($N) = @_; # name
	if ($self->{_sequences}{$N}) {
		return length($self->{_sequences}{$N});
	} else {return 0}
}

sub getNucleotide
{
	my $self = shift;
	my ($N,$P) = @_; # name, position
	return $self->getSequence($N, $P, $P);
}

sub getSequence
{
	my $self = shift;
	my ($N, $S, $E) = @_; # NSE identifier: Name, Start, End (1-based)
	my $L = $E ? $E-$S+1 : '';
	if ($self->{_sequences}{$N} && $L ne '') {
		return substr($self->{_sequences}{$N}, $S-1, $L);
	} elsif ($self->{_sequences}{$N} && $L eq '') {
		return substr($self->{_sequences}{$N}, $S-1);
	} else {return 0}
}

sub getSequenceSubstr
{
	my $self = shift;
	my ($N, $S0, $L) = @_; # Just like Perl substr function: 0-based start and length
	if ($self->{_sequences}{$N}) {
		return substr($self->{_sequences}{$N}, $S0, $L);
	} else {return 0}
}

sub setSequence
{
	my $self = shift;
	my ($N, $S, $E, $insert) = @_; # Name, Start (1-based), End, sequence to insert
	my $L = $E-$S+1;
	substr($self->{_sequences}{$N}, $S-1, $L) = $insert;
}

sub findSequence
{
	my $self = shift;
	my ($seq) = @_;
	my $length = length $seq;
	my %results; # {name}{start} = end
	for my $N (keys %{$self->{_sequences}}) {
		my $offset = index($self->{_sequences}{$N}, $seq, 0);
		if ($offset == -1) {next}
		my $end = $offset + $length;
		$offset++;
		$results{$N}{$offset} = $end;
		while ($offset != 0) {
			$offset = index($self->{_sequences}{$N}, $seq, $offset);
			$end = $offset + $length;
			$offset++;
			unless ($offset == 0) {$results{$N}{$offset} = $end}
		}
	}
	if (!%results) {return 0}
	else {return \%results}
}

sub freeSequence
{
	my $self = shift;
	my ($N) = @_;
	$self->{_sequences}{$N} = '';
	return $N;
}

## PRIVATE METHODS ##
sub _load
{
	my $self = shift;
	my $current_header = 0; # set to false in order to test for proper header assignment before loading sequence
	if (-d$self->{_fasta_location}) {
		if ($self->{_messages} > 0) {print STDERR "Loading FASTAs...\n"}
		my $dir = $self->{_fasta_location};
		opendir(FASTADIR, $dir);
		my @dir_files = readdir(FASTADIR);
		closedir(FASTADIR);
		for my $filename (@dir_files) {
			if ($filename=~/\.fa$/i || $filename=~/\.fasta$/i) {
				open FASTA, "<$dir/$filename" || die "[QueryFASTA] Could not open FASTA: $!\n";
				$self->{_fasta_count}++;
				while (<FASTA>) {
					if (/^\#/ || /^\;/) {next}
					chomp;
					if (/^>/) {
						($current_header = $_) =~ s/^>//;
						$current_header = (split(/\s/, $current_header))[0];
						$current_header =~ s/\|.+$//;
						if ($self->{_messages} > 1) {print STDERR "$current_header\n"}
					} else {
						if ($current_header) {
							$self->{_sequences}{$current_header} .= $_;
						}
					}	
				}
				close FASTA;
			}
		}
	} else {
		if ($self->{_messages} > 0) {print STDERR "Loading FASTAs...\n"}
		my $filename = $self->{_fasta_location};
		open FASTA, "<$filename" || die "[QueryFASTA] Could not open FASTA: $!\n";
		$self->{_fasta_count}++;
		while (<FASTA>) {
			if (/^\#/ || /^\;/) {next}
			chomp;
			if (/^>/) {	
				($current_header = $_) =~ s/^>//;
				$current_header = (split(/\s/, $current_header))[0];
				$current_header =~ s/\|.+$//;
				if ($self->{_messages} > 1) {print STDERR "$current_header\n"}
			} else {
				if ($current_header) {
					$self->{_sequences}{$current_header} .= $_;
				}
			}	
		}
		close FASTA;
	}
	return $self->{_fasta_count};
}

1;
