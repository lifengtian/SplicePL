#!/usr/bin/perl


### Find unique junctions
### Lifeng Tian
### Apr 2010

my $usage = "perl unique_junctions.pl [junction file]";
while(<>){
if (/^chr/ ){
	my @a=split/\t/;
	$h{$a[0]}++;
	$line{$a[0]}=$_;
}
}
foreach (sort keys %h){
		print $line{$_};
}
