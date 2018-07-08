use strict ;
use warnings ;

my$in = $ARGV[0] ;
my$out1 = $in ;
my$out2 = $in ;
$out1 =~ s/(.+)_coverage_by_site.txt/$1_aberrant_regions.txt/ ;
$out2 =~ s/(.+)_coverage_by_site.txt/$1_confident_regions.txt/ ;

my$scaff ;
my$site ;
my$depth ;

my%coverage ;
my%total ;

my%window1 ;
my%window1_count  ;
my%num_window1 ;

my%window2 ;
my%window2_count ;
my%num_window2 ;

open IN, "<$in" ;

while (<IN>) {

    if ( $_ =~ m/^#/ ) {
        next ;
    }

    my@split = split(/\t/, $_) ;

    $scaff = $split[0] ;
    $site = $split[1] ;
    $depth = $split[2] ;

    if (!$window1_count{$scaff}) {
        $window1_count{$scaff} = 0 ;
        $window2_count{$scaff} = -500 ;
        }

    if (!$num_window1{$scaff}) {
        $num_window1{$scaff} = 1 ;
        $num_window2{$scaff} = 1 ;
        }

    $coverage{$scaff} += $depth ;
    $total{$scaff} ++ ;

    if ($window1_count{$scaff} < 1000) {
        $window1{$scaff}{$num_window1{$scaff}}{"TOTAL_DEPTH"} += $depth ;
        push @{$window1{$scaff}{$num_window1{$scaff}}{"SITES"}}, $site ;
        $window1_count{$scaff} ++ ;
        }

    else {
        $num_window1{$scaff} ++ ;

        $window1{$scaff}{$num_window1{$scaff}}{"TOTAL_DEPTH"} += $depth ;
        push @{$window1{$scaff}{$num_window1{$scaff}}{"SITES"}}, $site ;
        $window1_count{$scaff} = 1 ;
        }

    if ($window2_count{$scaff} < 0 ) {
        $window2_count{$scaff} ++ ;
        }

    elsif ($window2_count{$scaff} >=0 && $window2_count{$scaff} < 1000) {
        $window2{$scaff}{$num_window2{$scaff}}{"TOTAL_DEPTH"} += $depth ;
        push @{$window2{$scaff}{$num_window2{$scaff}}{"SITES"}}, $site ;
        $window2_count{$scaff} ++ ;
        }

    else {
        $num_window2{$scaff} ++ ;

        $window2{$scaff}{$num_window2{$scaff}}{"TOTAL_DEPTH"} += $depth ;
        push @{$window2{$scaff}{$num_window2{$scaff}}{"SITES"}}, $site ;
        $window2_count{$scaff} = 1 ;
        }
    }

close IN ;

my%aberrant ;
my%confident ;
my$avg_cov ;


foreach my$scaff (keys %coverage) {
    $avg_cov = $coverage{$scaff} / $total{$scaff} ;

    foreach my$window (keys %{$window1{$scaff}}) {
        my@sites = @{$window1{$scaff}{$window}{"SITES"}} ;

        my$window_avg = $window1{$scaff}{$window}{"TOTAL_DEPTH"} / $#sites ;

        if ($window_avg >= 2*$avg_cov || $window_avg <= $avg_cov / 2) {
            $aberrant{$scaff}{$window.".1"}{"START"} = $sites[0] ;
            $aberrant{$scaff}{$window.".1"}{"END"} = $sites[-1] ;
            $aberrant{$scaff}{$window.".1"}{"COV"} = $window_avg ;
        }

        else {
            $confident{$scaff}{$window.".1"}{"START"} = $sites[0] ;
            $confident{$scaff}{$window.".1"}{"END"} = $sites[-1] ;
            $confident{$scaff}{$window.".1"}{"COV"} = $window_avg ;
        }
    }

    foreach my$window (keys %{$window2{$scaff}}) {
        my@sites = @{$window2{$scaff}{$window}{"SITES"}} ;

        my$window_avg = $window2{$scaff}{$window}{"TOTAL_DEPTH"} / $#sites ;

        if ($window_avg >= 2*$avg_cov || $window_avg <= $avg_cov / 2) {
            $aberrant{$scaff}{$window.".2"}{"START"} = $sites[0] ;
            $aberrant{$scaff}{$window.".2"}{"END"} = $sites[-1] ;
            $aberrant{$scaff}{$window.".2"}{"COV"} = $window_avg ;
        }

        else {
            $confident{$scaff}{$window.".2"}{"START"} = $sites[0] ;
            $confident{$scaff}{$window.".2"}{"END"} = $sites[-1] ;
            $confident{$scaff}{$window.".2"}{"COV"} = $window_avg ;
        }
    }
}

open OUT1, ">$out1" ;

foreach my$scaff (sort keys %aberrant) {
    foreach my$window (sort {$a<=>$b} keys %{$aberrant{$scaff}}) {
        print OUT1 $scaff, "\t", $avg_cov, "\t", $window, "\t", $aberrant{$scaff}{$window}{"START"}, "\t", $aberrant{$scaff}{$window}{"END"}, "\t", $aberrant{$scaff}{$window}{"COV"}, "\n" ;
        }
    }

close OUT1 ;


open OUT2, ">$out2" ;

foreach my$scaff (sort keys %confident) {
    foreach my$window (sort {$a<=>$b} keys %{$confident{$scaff}}) {
        print OUT2 $scaff, "\t", $avg_cov, "\t", $window, "\t", $confident{$scaff}{$window}{"START"}, "\t", $confident{$scaff}{$window}{"END"}, "\t", $confident{$scaff}{$window}{"COV"}, "\n" ;
        }
    }

close OUT2 ;