# SMD workflow (Russell-et-al-2018)

Scripts used in Russell S, McCartney E, and Cavanaugh C. Evidence of transovarial symbiont transmission in the chemosynthetic bivalve Solemya velum. In prep.


How to use this workflow to identify candidate genes for detection of a specific microbial group via PCR/qPCR/FISH:

1) For your taxon of interest, obtain the genome sequence in fasta format and the annotation in gff format.

2) Run extract_CDS.pl (in ./scripts) to extract the nucleic acid sequence for each of the coding regions annotated in the gff file:
   Usage: perl extract_cds.pl genome.fasta genome_cds.gff

3) Blast each gene in the genome against the NCBI non-redundant nucleotide database (see example command line blast script: ./test_example/Svelum_sym_genes-blastn, where the input file is the output file from step 2, specifically a fasta containing separate nucleic acid sequences for each gene)

4) Find candidate marker genes from the blast output by running marker_candidates.pl (in ./scripts):
  Usage: perl marker_candidates.pl blastn_hits.txt  
  

Optional: If population-level whole-genome data is available, the following steps can be run to produce a filtered candidate list that takes into account genome structure uncertainty (by filtering out regions with coverage anomolies) and population-level polymorphism that could hinder the design of primers that work across genotypes within the taxon of interest 

5) Run find_aberrant_regions.pl (in ./scripts) on each sample from the population to find regions of the genome with normal/abnormal coverage relative to the genome-wide average. The "coverage_by_site.txt" file is a tab-delimited file with one row for each position in the genome, with the first column containing the scaffold, the second column containing the position in bp, and the third column containing the mapped read coverage. A file of this type can be produced with Bedtools using the following command:
  Usage: genomeCoverageBed -d -ibam sample_realigned.bam -g genome.fasta > sample_coverage_by_site.txt

Then find the regions that match and mismatch the genome-wide average
  Usage: perl find_aberrant_regions.pl sample_coverage_by_site.txt

Put all the output "_confident_regions.txt" files in one directory, update the following script (in ./scripts) with your genome's scaffold lengths, and run the script to output the regions of the genome that were found to be of average coverage in all samples (output file == confident_ranges.txt):
  Usage: perl consensus_confident.pl ./confident_regions/

6) Use a variant calling algorithm such as GATK's UnifiedGenotyper to call variants in each sample relative to the reference genome sequence, and output a .vcf file containing variant sites for all samples

7) Filter candidate marker genes based on conserved regions of the genome (in regards to coverage and polymorphism) by running filter_candidates.pl (in ./scripts):
  perl filter_candidates.pl confident_ranges.txt Svelum_sym_genes_cds_candidate_markers.txt genome_cds.gff genome_variants.vcf
  

Next steps: Primers can be designed for the genes listed in "candidate_marker.txt" or "candidate_markers-filtered.txt", starting with the best candidates, at the top of the list, first. Gene lists are sorted in descending order by the number of hits reported for that gene, as this is a proxy for the data available on the gene's specificity. Candidate lists are tab delimted files, so additional sorting by percent identity, and number of variants (if following steps 5-7) can be performed in a text editor or spreadsheet.
