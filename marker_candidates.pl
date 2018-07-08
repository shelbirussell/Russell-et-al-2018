use strict ;
use warnings ;

my$blast = $ARGV[0] ;

my$sample = $blast ;
$sample =~ s/(.+)_.+/$1/ ;
my$output = $sample . "_candidate_markers.txt" ;

open BLAST, "<$blast" ;

my%hits ;
my%counts ;
my%pidents ;

while (<BLAST>) {
    if ($_ =~ m/^#/) {
        next ;
        }

    my@hit = split(/\t/, $_) ;

    my$qseqid = $hit[0] ;
    my$sallacc = $hit[1] ;
    my$salltitles = $hit[2] ;
    my$pident = $hit[3] ;
    my$align_length = $hit[4] ;
    my$mismatch = $hit[5] ;
    my$gapopen = $hit[6] ;
    my$qstart = $hit[7] ;
    my$qend = $hit[8] ;
    my$sstart = $hit[9] ;
    my$send = $hit[10] ;
    my$evalue = $hit[11] ;
    my$bitscore = $hit[12] ;

    $hits{$qseqid}{$sallacc}{"TITLES"} = $salltitles ;
    $hits{$qseqid}{$sallacc}{"PIDENTITY"} = $pident ;
    $hits{$qseqid}{$sallacc}{"LENGTH"} = $align_length ;
    $hits{$qseqid}{$sallacc}{"MISMATCH"} = $mismatch ;
    $hits{$qseqid}{$sallacc}{"GAPS"} = $gapopen ;
    $hits{$qseqid}{$sallacc}{"QSTART"} = $qstart ;
    $hits{$qseqid}{$sallacc}{"QEND"} = $qend ;
    $hits{$qseqid}{$sallacc}{"SSTART"} = $sstart ;
    $hits{$qseqid}{$sallacc}{"SEND"} = $send ;
    $hits{$qseqid}{$sallacc}{"Evalue"} = $evalue ;
    $hits{$qseqid}{$sallacc}{"BITSCORE"} = $bitscore ;

    if (!$counts{$qseqid}) {
        $counts{$qseqid} = 1 ;
        }

    else {
        $counts{$qseqid} ++ ;
        }

    $pidents{$qseqid}{"TOTAL"} += $pident ;

}

close BLAST ;

my%avg_pidents ;

foreach my$query (keys %pidents) {
    $avg_pidents{$query} = $pidents{$query}{"TOTAL"} / $counts{$query} ;
    }

my@rank_hit_abundance = sort { $counts{$b} <=> $counts{$a} } keys %counts ;

open OUT, ">$output" ;

foreach my$hit (@rank_hit_abundance) {
    if ($counts{$hit} < 100) {
        next ;
        }

    if ($avg_pidents{$hit} > 80 ) {
        next ;
        }

    print OUT $hit, "\t", $counts{$hit}, "\t", $avg_pidents{$hit}, "\n" ;
    }