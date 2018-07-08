use strict ;
use warnings ;
use Getopt::Long ;

my$blast ;
my$gff ;
my$vcf ;
my$regions ;

GetOptions ('blast=s' => \$blast,
			'gff=s' => \$gff,
			'vcf=s' => \$vcf,
			'regions=s' => \$regions) ;
or die ("Error in command line arguments\n") ;

