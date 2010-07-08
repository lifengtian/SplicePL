#!/usr/bin/perl

use strict;
use FileHandle;

=head1 Note
The first three required BED fields are:
   1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random).
   2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
   3. chromEnd - The ending position of the feature in the chromosome or scaffold. 
   The chromEnd base is NOT included in the display of the feature. 
   For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 


The first column in the junction file is formatted as chrx:start-end, e.g., chr1:1-100
    Here, the first base in a chromosome is numbered 1!
    
=cut

my $usage = "perl junction_to_bed [junction file] [bed file]";

if ( @ARGV != 2 ) {
		die $usage;
}

my ($junctionFileName, $bedFileName) = @ARGV;
my $junctionFileHandle = FileHandle->new($junctionFileName,"r");
my $bedFileHandle = FileHandle->new($bedFileName,"w");

while ( <$junctionFileHandle> ){
	my $j = (split/\t/)[0];
	my ($chr, $start,$end) = split/[:-]/,$j;
	print $bedFileHandle join("\t",$chr, $start - 1, $end), "\n";
}

undef $junctionFileHandle;
undef $bedFileHandle;