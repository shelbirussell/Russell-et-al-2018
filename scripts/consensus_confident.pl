use strict ;
use warnings ;
use List::Compare ;

## This script finds the regions confident in all samples with coverage above (10x).
## Deactivated: If region is not supported by only KB25, the region is allowed to pass to the consensus because KB25 has apparent coverage problems.
## Usage: perl consensus_confident.pl ./confident_regions/

my$dir = $ARGV[0] ;

## This is specific to the S. velum symbiont/mito genomes
## Replace with the lengths of the contigs in your species' genome fasta file
my %ends ;
$ends{"Sv_mito_chromosome"} = 15660 ;
$ends{"Sv_sym_scaffold1"} = 1213831 ;
$ends{"Sv_sym_scaffold2"} = 892555 ;
$ends{"Sv_sym_scaffold3"} = 537613 ;
$ends{"Sv_sym_scaffold4"} = 28016 ;

opendir(SOURCE, $dir) ;
my@files = readdir(SOURCE) ;

my%regions ;

my@names ;
foreach my$file (@files) {
    if ($file =~ m/.+confident_regions.txt/) {
        $file =~ m/^(\w+\d+gill).+/ ;
        push @names, $1 ;
        }
    }

foreach my$file (@files) {
    if ($file =~ m/.+confident_regions.txt/) {

        open IN, "<$dir$file" ;

        my$sample = $file ;
        $sample =~ s/^(\w+\d+gill).+/$1/ ;

        while (<IN>) {

            if ($_ =~ m/^#/) {
                next ;
                }

            my@split = split(/\t/, $_) ;
            my$scaff = $split[0] ;
            my$avg_cov = $split[1] ;
            my$window = $split[2] ;

            if($avg_cov < 10){
                @names = grep {$_ ne $sample} @names ;
                last ;
                }

            push @{$regions{$scaff}{$window}}, $sample ;

            }

        close IN ;

    }

}

my%confident ;

foreach my$scaff (keys %regions) {
    foreach my$window (sort {$a<=>$b} keys %{$regions{$scaff}}) {
        my@samples = @{$regions{$scaff}{$window}} ;

        if ($#samples >= $#names) {
            push @{$confident{$scaff}}, $window ;
            }

        else {
            my$lc = List::Compare->new(\@samples, \@names) ;
            my@unique = $lc->get_symmetric_difference() ;

#            if (scalar @unique == 1 && "KB25gill" ~~ @unique) {
#                push @{$confident{$scaff}}, $window ;
#                }

#            else {
#                print $#samples, " ", $#names, "\t", $scaff, "\t", $window, "\t", "absent_samples:", join(",", @unique), "\n";
#                next ;
#                }
            }
        }
    }

my%ranges ;

foreach my$scaff (keys %confident) {
    my@windows = @{$confident{$scaff}} ;

    foreach my$window (@windows) {
        my$range_end ;
        my$range_start ;

        if ($window =~ m/\d+\.1/) {
            my$win_count = $window ;
            $win_count =~ s/(\d+)\.1/$1/ ;

            if ($win_count * 1000 > $ends{$scaff}) {
                $range_end = $ends{$scaff} ;
                $range_start = $win_count * 1000 - 999 ;

                last ;
            }

            else {
                $range_end = $win_count * 1000 ;
                $range_start = $range_end - 999 ;
            }
        }

        if ($window =~ m/\d+\.2/) {
            my$win_count = $window ;
            $win_count =~ s/(\d+)\.2/$1/ ;

            if ($win_count * 1000 + 500 > $ends{$scaff}) {
                $range_end = $ends{$scaff} ;
                $range_start = $win_count * 1000 - 499 ;

                last ;
            }

            else {
                $range_end = $win_count * 1000 + 500 ;
                $range_start = $range_end - 999 ;
            }
        }

        if (!$ranges{$scaff}) {
            my@start_range = ($range_start .. $range_end) ;
            $ranges{$scaff} = \@start_range ;
        }

        else {
            my@add_range ;
            my@scaff_ranges = @{$ranges{$scaff}} ;

            if ($range_start ~~ @scaff_ranges) {

                if ($range_end ~~ @scaff_ranges) {
                    next ;
                }

                else {
                    @add_range = $scaff_ranges[-1]+1 .. $range_end ;
                    push @{$ranges{$scaff}}, @add_range ;
                }
            }

            else {
                @add_range = $range_start .. $range_end ;
                push @{$ranges{$scaff}}, @add_range ;
            }
        }
    }
}

foreach my$scaff (sort {$a cmp $b} keys %ranges) {
    my@range = @{$ranges{$scaff}} ;

    print $scaff, "\t", $range[0] ;

    foreach my$i (1..$#range) {
        my$previous = $i - 1 ;

        if ($range[$i] - 1 == $range[$previous]) {

            if ($range[$i] == $range[-1]) {
                print "-", $range[-1] ;
                }

            else {
                next ;
                }
            }

        else {
            print "-", $range[$previous], ",", $range[$i] ;
            }
        }

    print "\n" ;

    }
