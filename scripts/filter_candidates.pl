use strict ;
use warnings ;

## This script takes filters the marker candidates for detecting free-living Sv symbionts
## and removes genes not in confident regions of the reference genome (found with find_aberrant_regions.pl).
## Then, the number of variant sites within each gene are reported by referencing the haploid vcf file.
## Usage: perl filter_candidates.pl ../filter_ref/confident_ranges.txt Svelum_sym_genes_cds_candidate_markers.txt ~/Reference/Solemya_velum_MitoSym_top4.gff ../GATK/Sv_MS_stampy_UGploidy1_VSboth-filtered.filtered.vcf

my$regions = $ARGV[0] ;
my$candidates = $ARGV[1] ;
my$gff = $ARGV[2] ;
my$vcf = $ARGV[3] ;

open REGIONS, "<$regions" ;

my%confident ;

while (<REGIONS>) {

    if ($_ =~ m/^#/) {
        next ;
    }

    my@split = split(/\t/, $_) ;
    my$scaff = $split[0] ;
    my@hyphen_ranges = split(",", $split[1]) ;
    my@ranges ;

    foreach my$range (@hyphen_ranges) {
        my@bound = split("-", $range) ;
        my@whole_range = $bound[0]..$bound[1] ;
        push @ranges, @whole_range ;
    }

    $confident{$scaff} = \@ranges ;
}

close REGIONS ;

open LOCI, "<$candidates" ;

my%candidates ;

while (<LOCI>) {

    if ($_ =~ m/^#/) {
        next ;
    }

    chomp $_ ;

    my@split = split(/\t/, $_) ;
    my$gene = $split[0] ;
    my$hits = $split[1] ;
    my$perc_id = $split[2] ;

    $candidates{$gene}{"HITS"} = $hits ;
    $candidates{$gene}{"IDENTITY"} = $perc_id ;

}

close LOCI ;

open GFF, "<$gff" ;

while (<GFF>) {

    if ($_ =~ m/^#/) {
        next ;
    }

    chomp $_ ;
    my @split = split ( /\t/, $_ ) ;

    if ( $split[2] ne "CDS" ) { next ; }

    my $gff_gene = "" ;
    if ( $split[8] =~ m/locus_tag=(.+)(;product|;note)/ ) {
        $gff_gene = $1 ;
    }
    if ( $split[8] =~ m/locus_tag=(.+);note/ ) {
        $gff_gene = $1 ;
    }

    if (exists $candidates{$gff_gene}) {
        $candidates{$gff_gene}{"SCAFF"} = $split[0] ;
        $candidates{$gff_gene}{"START"} = $split[3] ;
        $candidates{$gff_gene}{"END"} = $split[4] ;
    }
}

close GFF ;


open VCF, "<$vcf" ;

my%variant_sites ;

while (<VCF>) {
    if ($_ =~ m/^##/) { next ; }

    my @split = split(/\t/, $_) ;

    my$scaff = $split[0] ;
    my$pos = $split[1] ;
    my$ref = $split[3] ;
    my$alt = $split[4] ;

    push @{$variant_sites{$scaff}}, $pos ;

    }

close VCF ;

my%filtered_candidates ;
my%hits ;

foreach my$gene (keys %candidates) {
    my$scaff = $candidates{$gene}{"SCAFF"} ;

    if ($candidates{$gene}{"START"} ~~ @{$confident{$scaff}} && $candidates{$gene}{"END"} ~~ @{$confident{$scaff}}) {
        $filtered_candidates{$gene}{"HITS"} = $candidates{$gene}{"HITS"} ;
        $filtered_candidates{$gene}{"IDENTITY"} = $candidates{$gene}{"IDENTITY"} ;
        $filtered_candidates{$gene}{"SCAFF"} = $candidates{$gene}{"SCAFF"} ;
        $filtered_candidates{$gene}{"START"} = $candidates{$gene}{"START"} ;
        $filtered_candidates{$gene}{"END"} = $candidates{$gene}{"END"} ;

        $hits{$gene} = $candidates{$gene}{"HITS"} ;
    }
}

my@rank_hit_abundance = sort { $hits{$b} <=> $hits{$a} } keys %hits ;

print "#Gene\tnumber_hits\tpercent_identity\tnumber_variants\n" ;

foreach my$gene (@rank_hit_abundance) {
    my$scaff = $filtered_candidates{$gene}{"SCAFF"} ;
    my$start = $filtered_candidates{$gene}{"START"} ;
    my$end = $filtered_candidates{$gene}{"END"} ;
    my@range = $start..$end ;

    my$variant_count = 0 ;

    foreach my$position (@range) {
        if ($position ~~ @{$variant_sites{$scaff}}) {
            $variant_count ++ ;
            }
        }

    print $gene, "\t", $filtered_candidates{$gene}{"HITS"}, "\t", $filtered_candidates{$gene}{"IDENTITY"}, "\t", $variant_count, "\n" ;

    }

