#!/usr/bin/perl

use strict;
### compare two  splice junction files
### A junctions file has its first column like chr1:1-20
### INPUT: two junction files
### OUTPUT: stdout each line will have 
###         1. ONE     only in first input file
###         2. TWO     only in second input file
###         3. BOTH    in both first and second input file
###



my $usage="usage: perl compare_two_junction_files.pl <junction1> <junction2>";
my $in1=$ARGV[0] || die $usage;
my $in2=$ARGV[1] || die $usage; 

open my $J1, $in1 or die "Error open $in1";
open my $J2, $in2 or die "Error open $in2";

my $count = 0;
my %junction1Lines;
my %junctionCounter;
while(<$J1>){
    chomp;
	my @a=split/\t/;
	if ($a[0]=~/^chr.+:\d+-\d+/){
	$junction1Lines{$a[0]}=$_;
	$junctionCounter{$a[0]}++;
	} else {
	die  "Exit! Wrong junction format in line $. ", $_,"\n";
	}
}

my %junction2Lines;
while (<$J2>){
    chomp;
    my @a=split/\t/;
        if ($a[0]=~/^chr.+:\d+-\d+/){
        $junction2Lines{$a[0]}=$_;
        $junctionCounter{$a[0]}++;
        } else {
        die   "Exit! Wrong junction format in line $. ", $_,"\n";
        }

}


foreach (sort keys %junctionCounter){
    if ( $junctionCounter{$_} == 2 ) {
	print $junction1Lines{$_},"\tBOTH\n";
    } elsif ( $junctionCounter{$_} == 1 ) {
	if ( $junction1Lines{$_} ) {
	    print $junction1Lines{$_},"\tONE\n";
	} else {
	    print $junction2Lines{$_},"\tTWO\n";
	}

    } else {
	print STDERR $junctionCounter{$_},"\tunknown\n";
	}
}

undef $J1;
undef $J2;

