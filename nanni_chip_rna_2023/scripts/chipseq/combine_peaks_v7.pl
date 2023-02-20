#! /usr/bin/perl


# peaks called from individual replicates were stacked
# "F" fusions with peaks from all replicates will be retained


use strict;
use warnings;
use POSIX;
use Getopt::Long;

my $max = '';
my $min = '';
my $m2m = '';
my $acf = '';
my $acaw = '';
my $summit = '';
my $asf = '';
my $asaw = '';
my $num_reps = '';
my $name_count = 1;
my $files = '';
my %storage;

GetOptions('max' => \$max, 'min'=> \$min, 'm2m' => \$m2m, 'acf=i' => \$acf, 'acaw' => \$acaw, 'summit' => \$summit, 'asf=i' => \$asf, 'asaw' => \$asaw, 'reps=i' => \$num_reps, 'files=s' => \$files);

my @file_list = split(/,/,$files);

my %peak_list;
my %rep_list;

foreach(@file_list){
    open(IN,"<","$_");
    my $sample = $_; $sample =~ s/\.sbed//g;

    $rep_list{$sample} = 0;
    


    while (<IN>) {
        chomp;
        my ($chr, $start, $end, $sumit) = split (/\t/, $_);
        my $line=$.-1;
        my $id="peak$line";
        my $midpoint = &calc_midpoint($start,$end);

        push( @{$peak_list{$chr}},{'id' => $id, 'start' => $start, 'end' => $end, 'sample' => $sample, 'mid' => $midpoint, 'summit' => $sumit} );

    }
    close(IN);
}

foreach my $chr (keys %peak_list){
    my @sorted_end =  sort { $a->{end} <=> $b->{end} } @{$peak_list{$chr}} ;
    my @sorted_start =  sort { $a->{start} <=> $b->{start} } @sorted_end ;

    if($max){
        my $hash = &find_max_overlap($chr,$num_reps,@sorted_start);
        %storage = (%{$hash},%storage);
    }
    if($min){
        my $hash = &find_min_overlap($chr,$num_reps,@sorted_start);
        %storage = (%{$hash},%storage);
    }
    if($m2m){
        my $hash = &find_m2m_overlap($chr,$num_reps,@sorted_start);
        %storage = (%{$hash},%storage);
    }
    if($acf){
        my $hash = &find_acf_overlap($chr,$num_reps,$acf,@sorted_start,);
        %storage = (%{$hash},%storage);
    }
    if($acaw){
        my $hash = &find_acaw_overlap($chr,$num_reps,@sorted_start,);
        %storage = (%{$hash},%storage);
    }

    if($summit){
        my $hash = &find_summit_overlap($chr,$num_reps,@sorted_start);
        %storage = (%{$hash},%storage);
    }
    if($asf){
        my $hash = &find_asf_overlap($chr,$num_reps,$asf,@sorted_start,);
        %storage = (%{$hash},%storage);
    }
    if($asaw){
        my $hash = &find_asaw_overlap($chr,$num_reps,@sorted_start,);
        %storage = (%{$hash},%storage);
    }
}

my $rep_names;
for my $key (sort keys %rep_list){
    $rep_names .= ",$key";
}
print "peak_id,chrom,peak_start,peak_end$rep_names\n";

&print_peak_region(\%storage);
    
#===  FUNCTION  ================================================================
#         NAME: find_max_overlap
#      PURPOSE: This function will fuse peaks together and create an area of
#               maximum overlapp.
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the maximum peak
#      region
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_max_overlap {
    my ($chr,$reps,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{start};
            $next_end = $array[$j]{end};

            until ($next_start > $peak_end || $j == @array ) {

                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            $hash{$name_count}{'chr'}   = $chr;
            $hash{$name_count}{'start'} = $peak_start;
            $hash{$name_count}{'end'}   = $peak_end;
            $hash{$name_count}{'reps'}   = $rep_cat;
            $name_count++;
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_min_overlap
#      PURPOSE: This function will find the minumum overlapping region.
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               @array = sorted array of information from a BED file
#
#      RETURNS: Prints to STDOUT a csv in BED format
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================
sub find_min_overlap {
    my ($chr,$reps,@array) = @_;
    my $prev_start = 999999999999999999999999;
    my $prev_end = 999999999999999999999999;
    my %hash;
    my $count;
    my $i = 0;
    my $j = 0;

    while($i < @array){

        $count = 1;
        my $j = $i + 1;
        my $sum_reps = 0;

        my $max_start;
        my $min_end;
        my @ends = ();
        my @starts = ();
        my $next_start;
        my $next_end;

        $rep_list{$array[$i]{'sample'}} = 1;

        my $cur_start = $array[$i]{start};
        my $cur_end = $array[$i]{end};
        $max_start = $cur_start;
        $min_end = $cur_end;

        push(@starts,$cur_start);
        push(@ends,$cur_end);

        if($j < @array){
            $next_start = $array[$j]{start};
            $next_end = $array[$j]{end};

            until ($next_start > $min_end || $j == @array){

                $rep_list{$array[$j]{'sample'}} = 1;
                push(@starts,$next_start);
                push(@ends,$next_end);

                my @sort_starts = sort { $b <=> $a } @starts;
                $max_start = $sort_starts[0];

                my @sort_ends = sort { $a <=> $b } @ends;
                $min_end = $sort_ends[0];

                $j++;
                $count++;

                if($j < @array){
                    $next_start = $array[$j]{start};
                    $next_end = $array[$j]{end};
                }
            }

        }
        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($max_start != $prev_start && $min_end != $prev_end && $sum_reps >= $reps ){
            $hash{$name_count}{'chr'}   = $chr;
            $hash{$name_count}{'start'} = $max_start;
            $hash{$name_count}{'end'}   = $min_end;
            $hash{$name_count}{'reps'}   = $rep_cat;

            $name_count++;
            $prev_start = $max_start;
            $prev_end = $min_end;
        }

        $i++;
    }

    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_m2m_overlap
#      PURPOSE: This function will fuse overlapping peaks and identify the peak
#      region as the middle to middle distance between the two most extreme
#      peaks.
#
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the
#      middle-to-middle peak region.
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_m2m_overlap {
    my ($chr,$reps,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my @mid_array = $array[$i]{'mid'};
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{'start'};
            $next_end = $array[$j]{'end'};

            until ($next_start > $peak_end || $j == @array ) {

                push(@mid_array,$array[$j]{'mid'});
                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            if($count == 1){
                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $peak_start;
                $hash{$name_count}{'end'}   = $peak_end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
            else {
                my @sorted = sort {$a <=> $b} @mid_array;
                my $min_mid = shift @sorted;
                my $max_mid = pop @sorted;

                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $min_mid;
                $hash{$name_count}{'end'}   = $max_mid;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_summit2summit_overlap
#      PURPOSE: This function will fuse overlapping peaks and identify the peak
#      region as the summit to summit distance between the two most extreme
#      peaks.
#
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the
#      summit-to-summit peak region.
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_summit_overlap {
    my ($chr,$reps,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my @summit_array = $array[$i]{'summit'};
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{'start'};
            $next_end = $array[$j]{'end'};

            until ($next_start > $peak_end || $j == @array ) {

                push(@summit_array,$array[$j]{'summit'});
                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            if($count == 1){
                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $peak_start;
                $hash{$name_count}{'end'}   = $peak_end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
            else {
                my @sorted = sort {$a <=> $b} @summit_array;
                my $min_summit = shift @sorted;
                my $max_summit = pop @sorted;
                
                if($min_summit == $max_summit){
                    $min_summit-- ;
                }

                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $min_summit;
                $hash{$name_count}{'end'}   = $max_summit;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_acf_overlap
#      PURPOSE: This functions identifies the average center location and takes
#      a footprint value to create peak region.
#
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               $fp = footprint
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the
#      acf peak region.
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_acf_overlap {
    my ($chr,$reps,$acf,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my @mid_array = $array[$i]{'mid'};
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{'start'};
            $next_end = $array[$j]{'end'};

            until ($next_start > $peak_end || $j == @array ) {

                push(@mid_array,$array[$j]{'mid'});
                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            if($count == 1){
                my $mid = &calc_midpoint($peak_start,$peak_end);
                my $foot = POSIX::floor($acf/2);
                my $start = $mid - $foot;
                my $end = $mid + $foot;
                
                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
            else {
                my $mean_mid = POSIX::ceil(&calc_array_mean(@mid_array));
                my $foot = POSIX::floor($acf/2);
                my $start = $mean_mid - $foot;
                my $end = $mean_mid + $foot;

                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_asf_overlap
#      PURPOSE: This functions identifies the average summit location and takes
#      a footprint value to create peak region.
#
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               $fp = footprint
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the
#      asf peak region.
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_asf_overlap {
    my ($chr,$reps,$asf,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my @summit_array = $array[$i]{'summit'};
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{'start'};
            $next_end = $array[$j]{'end'};

            until ($next_start > $peak_end || $j == @array ) {

                push(@summit_array,$array[$j]{'summit'});
                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            if($count == 1){
                my $summit = $summit_array[0];
                my $foot = POSIX::floor($asf/2);
                my $start = $summit - $foot;
                my $end = $summit + $foot;
                
                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
            else {
                my $mean_summit = POSIX::ceil(&calc_array_mean(@summit_array));
                my $foot = POSIX::floor($asf/2);
                my $start = $mean_summit - $foot;
                my $end = $mean_summit + $foot;

                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_acaw_overlap
#      PURPOSE: This functions identifies the average center location and the
#      average peak width. It then uses these to measures to create the peak
#      region.
#
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the
#      acaw peak region.
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_acaw_overlap {
    my ($chr,$reps,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my @mid_array = $array[$i]{'mid'};
        my @foot_array = $array[$i]{'end'} - $array[$i]{'start'};
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{'start'};
            $next_end = $array[$j]{'end'};

            until ($next_start > $peak_end || $j == @array ) {

                my $width = $next_end - $next_start;
                push(@foot_array,$width);

                push(@mid_array,$array[$j]{'mid'});
                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            if($count == 1){
                my $mid = &calc_midpoint($peak_start,$peak_end);
                my $width = pop(@foot_array);
                my $foot = POSIX::floor($width/2);
                my $start = $mid - $foot;
                my $end = $mid + $foot;
                
                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
            else {
                my $mean_mid = POSIX::ceil(&calc_array_mean(@mid_array));
                my $mean_width = POSIX::ceil(&calc_array_mean(@foot_array));
                my $foot = POSIX::floor($mean_width/2);
                my $start = $mean_mid - $foot;
                my $end = $mean_mid + $foot;

                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: find_asaw_overlap
#      PURPOSE: This functions identifies the average summit location and the
#      average peak width. It then uses these to measures to create the peak
#      region.
#
#   PARAMETERS: $chr = chromosome
#               $reps = replicate number
#               @array = sorted array of information from a BED file
#
#      RETURNS: A hash containing peak_ID, chr, start, end for the
#      acaw peak region.
#
#     COMMENTS: It is important that the array is sorted first by end then by
#     start
#
#===============================================================================

sub find_asaw_overlap {
    my ($chr,$reps,@array) = @_;
    my %hash;

    for (my $i=0; $i < @array; $i++){
        my $peak_start = $array[$i]{start};
        my $peak_end = $array[$i]{end};
        my $next_start = 9999999999999999999999;
        my $next_end = 999999999999999999999;
        my $j = $i + 1;
        my $count = 1;
        my @summit_array = $array[$i]{'summit'};
        my @foot_array = $array[$i]{'end'} - $array[$i]{'start'};
        my $sum_reps = 0;
        $rep_list{$array[$i]{'sample'}} = 1;

        if($j < @array){
            $next_start = $array[$j]{'start'};
            $next_end = $array[$j]{'end'};

            until ($next_start > $peak_end || $j == @array ) {

                my $width = $next_end - $next_start;
                push(@foot_array,$width);

                push(@summit_array,$array[$j]{'summit'});
                my @ends = ($peak_end, $next_end);
                my @sorted = sort { $a <=> $b } @ends;

                $peak_end = pop @sorted;
                $rep_list{$array[$j]{'sample'}} = 1;

                $i++;
                $j = $i + 1;
                $count++;

                if ($j < @array){
                    $next_end = $array[$j]{end};
                    $next_start = $array[$j]{start};
                }
            } 
        }

        my $rep_cat = ''; 

        for my $key (sort keys %rep_list){
            $sum_reps += $rep_list{$key};
            $rep_cat .= ",$rep_list{$key}";
            $rep_list{$key} = 0;
        }

        if($sum_reps >= $reps){
            if($count == 1){
                my $summit = $summit_array[0];
                my $width = pop(@foot_array);
                my $foot = POSIX::floor($width/2);
                my $start = $summit - $foot;
                my $end = $summit + $foot;
                
                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
            else {
                my $mean_summit = POSIX::ceil(&calc_array_mean(@summit_array));
                my $mean_width = POSIX::ceil(&calc_array_mean(@foot_array));
                my $foot = POSIX::floor($mean_width/2);
                my $start = $mean_summit - $foot;
                my $end = $mean_summit + $foot;

                $hash{$name_count}{'chr'}   = $chr;
                $hash{$name_count}{'start'} = $start;
                $hash{$name_count}{'end'}   = $end;
                $hash{$name_count}{'reps'}   = $rep_cat;
                $name_count++;
            }
        }
    }
    return(\%hash);
}

#===  FUNCTION  ================================================================
#         NAME: print_peak_region
#      PURPOSE: Prints to STDOUT the peak region.
#   PARAMETERS: A Storage hash passed by reference containing all of the peak
#   values
#
#     COMMENTS: none
#
#===============================================================================
sub print_peak_region {
    my $storage = shift;
    my %storage_hash = %{$storage};
    
    for my $key (sort keys %storage_hash){
        print "peak$key,$storage_hash{$key}{'chr'},$storage_hash{$key}{'start'},$storage_hash{$key}{'end'}$storage_hash{$key}{'reps'}\n";
    }
}

#===  FUNCTION  ================================================================
#         NAME: calc_midpoint
#      PURPOSE: Calculates the midpoint of a start and end
#   PARAMETERS: $start = start point
#               $end = end point
#
#     COMMENTS: none
#
#===============================================================================
sub calc_midpoint {
    my ($start, $end) = @_;
    my $adj = POSIX::ceil(($end - $start)/2);
    my $mid = $start + $adj;
    return($mid);
}

#===  FUNCTION  ================================================================
#         NAME: calc_array_mean
#      PURPOSE: Calculates the mean of an array of number
#   PARAMETERS: @array = an array of values to take the mean.
#
#     COMMENTS: none
#
#===============================================================================
sub calc_array_mean {
    my (@array) = @_;
    my $sum = 0;
    $sum += $_ for(@array);
    my $mean = $sum / scalar(@array);
    return($mean);
}
