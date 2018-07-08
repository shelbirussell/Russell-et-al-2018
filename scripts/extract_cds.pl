use strict ;
use warnings ;
use File::Basename ;

## This script extracts gene sequences from a fasta using the coordinates in the associated gff file
## The genes are written to separate fasta files, with header information from the gff file
## usage: perl ../extract_cds.pl ~/Reference/Solemya_velum_MitoSym_top4.fasta ../Solemya_velum_MitoSym_top4confident_regions.gff

my$fasta = $ARGV[0] ;
my$gff = $ARGV[1] ;

my%seqs = read_fasta($fasta) ;
my%anno = read_gff($gff) ;

my$cds = 0 ;
foreach my$scaff (keys %anno) {
	foreach my$id (keys %{$anno{$scaff}}) {
		$cds ++ ;
	}
}

print "Number of coding sequences: ", $cds, "\n" ;

my$out = basename($fasta) ;
$out =~ s/.fasta// ;
open OUT, ">${out}_cds_nt.fasta" ;

foreach my$scaffold (keys %anno) {
    foreach my$gene (keys %{$anno{$scaffold}}) {

        my$length = abs($anno{$scaffold}{$gene}{"START"} - $anno{$scaffold}{$gene}{"STOP"}) + 1 ;
        my$cut = substr($seqs{$scaffold}, $anno{$scaffold}{$gene}{"START"}-1, $length) ;

        if ($anno{$scaffold}{$gene}{"STRAND"} eq "-" ) {reverse_complement($cut) ;}

        print OUT ">", $gene, "\t", $scaffold, ":", $anno{$scaffold}{$gene}{"START"}, "-", $anno{$scaffold}{$gene}{"STOP"}, "\tlength:", $length, "\tstrand:", $anno{$scaffold}{$gene}{"STRAND"}, "\tproduct:", $anno{$scaffold}{$gene}{"GENE"}, "\n" ;
        print OUT "$_\n" foreach ($cut =~ /.{1,80}/g) ;

    }
}

sub read_fasta {
    open FASTA, "<$fasta" ;

    my%seqs ;
    my$header ;
    my$seq ;

    while (<FASTA>) {

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
                $seqs{$header} = $seq ;
                }

            $header = $_ ;
            $header =~ s/^>// ;
            $header =~ s/\s+$// ;

            $seq = "" ;
            }

        else {
            $_ =~ s/\s+//g ;
            $seq .= $_ ;
            }
        }

    close FASTA ;

    if ($seq) {
        $seqs{$header} = $seq ;
        }

    return %seqs ;
    }

sub read_gff {
    open GFF, "<$gff" ;

    my %data ;

    while (<GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\s+/, $_ ) ;

        if ( $split[0] eq "Sv_mito_chromosome" ) { next ; }

        if ( $split[2] ne "CDS" ) { next ; }

        my $id = "" ;
        if ( $split[8] =~ m/locus_tag=(.+)(;product|;note)/ ) {
            $id = $1 ;
        }
        if ( $split[8] =~ m/locus_tag=(.+);note/ ) {
            $id = $1 ;
        }

        my $gene = "" ;
        if ( $split[8] =~ m/product=([\(\)a\/,a-zA-Z0-9_-]+)/ ) {
            $gene = $1 ;
        }

        $data{$split[0]}{$id}{"START"} = $split[3] ;
        $data{$split[0]}{$id}{"STOP"} = $split[4] ;
        $data{$split[0]}{$id}{"STRAND"} = $split[6] ;
        $data{$split[0]}{$id}{"GENE"} = $gene ;

       }
    close GFF ;
    return %data ;
    }

sub reverse_complement {
    my$dna = shift ;

    my$revcomp = reverse($dna) ;

    $revcomp =~ tr/ACGTacgt/TGCAtgca/ ;

    return $revcomp ;
}
