#!/usr/bin/perl
use strict;
use FileHandle;

### compare predicted splice junctions with RefSeq refFlat junctions
### RefSeq (hg19) splice junction file is in the SplicePL/data/RefSeq_hg19.junc
### output : two files named known.junc and novel.junc will be created in the current folder


my $usage = "compare_with_RefSeq_junctions.pl [predicted junction file] [RefSeq junction file] \n";

die $usage if ( @ARGV != 2 );

my ($predictedJunctionFileName, $refSeqJunctionFileName) = @ARGV ;
unless ( -f $refSeqJunctionFileName ){
	die $refSeqJunctionFileName, "does not exist. ", $usage;
}


my $pf=FileHandle->new($predictedJunctionFileName, "r");
my $rf=FileHandle->new($refSeqJunctionFileName,"r");
my $o_known=FileHandle->new("known.junc", "w");
my $o_novel=FileHandle->new("novel.junc","w");

my %ref_junc ;

# RefSeq junctions
while(<$rf>){
    my @a=split/\t/;
    $ref_junc{$a[0]}++;
}

# predicted junctions
my $numberOfRefSeqJunctions = 0;
my $numberOfNonRefSeqJunctions = 0;


while(<$pf>){
	chomp;
    my @a=split/\t/;
    if ( $ref_junc{$a[0]} ) {
		$numberOfRefSeqJunctions++;
		print $o_known $_,"\tREFSEQ\n";
    } else {
    	$numberOfNonRefSeqJunctions++;
		print $o_novel $_,"\tNOVEL\n";
    } 
}

print "RefSeq junctions: $numberOfRefSeqJunctions \n Novel junctions: $numberOfNonRefSeqJunctions \n";


undef $pf;
undef $rf;
undef $o_known;
undef $o_novel;


