#!/usr/bin/perl -w

# Perl module for SplicePL
#
# Please direct questions and support issues to <lifeng2@mail.med.upenn.edu> 
#
# Copyright Lifeng Tian
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

=head1 DESCRIPTION


=head1 Check List 
Before running the SplicePL, make sure:
1. BLAT binary is in the PATH environment variable.
2. SplicePL is in the PATH, and executable!
3. BLAT database file, i.e., the xxx.2bit file


=cut 


use strict;
use Carp;
use Cwd 'abs_path';
use Env '@PATH';

my $usage = "perl /path/to/SplicePL.pl <fasta file name> </path/to/database> <no. of threads> <tileSize> <stepSize> <minScore>\n";
if ( @ARGV != 6 ) {
	die $usage;
}







my $fasta = $ARGV[0];    # input fasta file
my $db = $ARGV[1];
my $chunk = $ARGV[2];    # number of threads to parallelize
my $tileSize = $ARGV[3]; # BLAT tileSize
my $stepSize = $ARGV[4]; # BLAT stepSize
my $minScore = $ARGV[5]; # BLAT minScore


############################### Check Executables, dababase etc. #############################

###### check the following binaries: BLAT

my $executable = 'blat';

unless ( grep -x "$_/$executable", @PATH ) {
     croak "Can't find blat executables in @PATH\n";
}

######  check the databases for BLAT
unless ( -f $db ) {
    croak "Can't find BLAT database file as $db\n";
}


######  current path #######
## You must run SplicePL in the same folder as fasta file
my $cp = abs_path();  

unless ( -f $cp.'/'.$fasta) {
    croak "Can't find $fasta file in $cp. Exit!";
}

if ( ! -d $cp.'/scripts' ) {
    mkdir( $cp.'/scripts');
}


################################# BLAT parameters #################################
my $repMatch = int(1024.0*$tileSize/$stepSize);
my $parameters="-noHead -repMatch=$repMatch -minScore=$minScore -stepSize=$stepSize -tileSize=$tileSize";
#my $parameters="-noHead -repMatch=$repMatch -minScore=$minScore -stepSize=$stepSize -tileSize=$tileSize -ooc=$pipeline/hg19.${tileSize}ooc.masked -out=pslx";




################################# Generate shell scripts to run Blat#################################
for ( my $i=1;$i<=$chunk;$i++) {
    my $command=<<EOF;
#!/bin/sh' 
date >> $cp/run.$i.log
echo blat $db $parameters $cp/$fasta.$i $cp/$fasta.$i.out >> $cp/run.$i.log
blat $db $parameters $cp/$fasta.$i $cp/$fasta.$i.out
date >> $cp/run.$i.log
EOF

	my $of='run.'.$i.'.sh';
	open OUT, ">$cp/scripts/".$of or die "Error open $of for writing. ";
	print OUT $command;
}

### Break up one input file into CHUNK pieces of input fasta file

&breakup_fasta($fasta,$chunk,$cp);



######  Confirm with user the pipeline status





### Run the BLAT pipeline
for my $p (1..$chunk) {
   my $pid = fork();
   if ($pid == -1) {
       die;
   } elsif ($pid == 0) {
    
    	my $blat_process="bash $cp/scripts/run.$p.sh ";
  	exec $blat_process or die;
   }
}
while (wait() != -1) {}
print "All $chunk BLAT jobs Done \n";




#################################### Run analysis #################################

## merge BLAT output files into one file
my $s='cat ';
my $temp = $cp.'/temp';
if ( ! -d $temp ) {
    mkdir($temp);
}

for my $j ( 1..$chunk) {
    $s .= $cp.'/'.$fasta.'.'.$j.'.out'.' ';
}

my $comm_cat=$s.'> '.$temp.'/'.$fasta.'.psl';
print $comm_cat,"\n";
system($comm_cat);



## mark reads status: U_pair, etc.
my $noMaxReads = &maxReads($cp.'/'.$fasta);

print  "noMaxReads ", $noMaxReads, "\n";

my $status = $temp.'/'.$fasta.'.psl'.'.status';
print $status,"\n";
&markReadsStatus ($noMaxReads, $temp.'/'.$fasta.'.psl', $status);

exit();


=head1

#u_pair
$find_status u_pair $status > $temp/u_pair
$get_u_pair_psl $temp/u_pair $temp/$fasta.psl > $temp/u_pair.psl
$check_str_u_pair $temp/u_pair.psl 1> $temp/u_pair_consistent.psl 2> $temp/u_pair_nonconsistent
date
echo Done 
#nu_pair
echo Find NU_PAIR status
date
$find_status nu_pair $status > $temp/nu_pair
$get_u_pair_psl  $temp/nu_pair $temp/$fasta.psl > $temp/nu_pair.psl
#$pslReps  -minNearTopSize=50 $temp/nu_pair.psl $temp/nu_pair_filtered.psl $temp/nu_pair_filtered.psr
cp $temp/nu_pair.psl $temp/nu_pair_filtered.psl
$tophit $temp/nu_pair_filtered.psl 1> $temp/nu_pair_filtered.tophit.psl 2> $temp/nu_pair_filtered.tophit.EXCLUDED.pslplusfirstcolumn
$mark_reads_status  $MAXREADS $temp/nu_pair_filtered.tophit.psl > $temp/nu_pair_tophit_status
$find_status u_pair $temp/nu_pair_tophit_status > $temp/nu_pair_tophit_u_pair

$get_u_pair_psl $temp/nu_pair_tophit_u_pair $temp/nu_pair_filtered.tophit.psl > $temp/u_pair_from_nu_pair.psl
$check_str_u_pair < $temp/u_pair_from_nu_pair.psl 1> $temp/u_pair_from_nu_pair_consistent.psl 2>$temp/u_pair_from_nu_pair_inconsistent_err
cat $temp/u_pair_consistent.psl $temp/u_pair_from_nu_pair_consistent.psl > $temp/u_pair_ALL.psl
date
echo Done
date
echo Done Prepare Qualified PSL output 


echo Find junctions
date

jp=$rp/junctions
mkdir $jp

$find_junctions $cp/$i/blat/temp/u_pair_ALL.psl > $jp/junctions
$uniquejunctions $jp/junctions > $jp/uniquejunctions
$filter_GTAG2 $jp/uniquejunctions | grep -h GOOD | awk '{if ($4>50) print $0}' > $jp/filtered_junctions
date
echo Done


=cut





#################################### Subroutines  #################################                                                                                                     

=head2
 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 
=cut



=head2 breakup_fasta
 Title   : breakup_fasta
 Usage   : &breakup_fasta($fasta, 6, "/path/to/blah")
 Function: split a fasta file into N pierces
 Returns : NULL. Side effect is number of chunks of files will be created
 Args    : Three scalars, first is the fasta file name,
           Second is the number of chunks
           Third is the current path
=cut

sub breakup_fasta {
    my ($fastafile, $numpieces,$cp) = @_;
    open(INFILE, $cp.'/'.$fastafile);

    print STDERR "counting fasta records in the $cp/$fastafile\n";

    my $filesize = `grep '>' $cp/$fastafile | wc -l`;
    chomp($filesize);
    $filesize =~ s/\s*//;
    print STDERR "No of fasta records: $filesize\n";

    my $numseqs = $filesize ;
    my $piecesize = int($numseqs / $numpieces);

    if ( $piecesize < 1 ) {
	die "Quit! Not enough pieces to break into $numpieces!";
    }    


    print STDERR "processing in $numpieces pieces of approx $piecesize size each.\n";
    if($piecesize % 2 == 1) {
        $piecesize++;
    }
    my $bflag = 0;
    for(my $i=1; $i<$numpieces; $i++) {
        my $outfilename = $cp.'/'.$fastafile . "." . $i;
        open(OUTFILE, ">$outfilename");
        for(my $j=0; $j<$piecesize; $j++) {
            my $line = <INFILE>;
            print OUTFILE $line;
	    $line=<INFILE>;
	    print OUTFILE $line;
        }
        close(OUTFILE);
    }
    my $outfilename = $cp.'/'.$fastafile. "." . $numpieces;
    open(OUTFILE, ">$outfilename");
    while(my $line2 = <INFILE>) {
        print OUTFILE $line2;
    }
    close(OUTFILE);
    return 0;
}





sub maxReads {
    my ($fn) = @_;
    my $r=`tail -n 2 $fn | head -1`;
    if ( $r=~/(\d+)/ ) {
	return $1;
    } else {
	croak "Did not find the maximal Reads number. Exit!";
    }
}




### mark_reads_status.pl <max readsID>
### Each read will be marked as
###    u_pair
###    u_forward
###    u_reverse
###    nu_pair
###    nu_forward
###    nu_reverse
###    nomap
###


### psl file has no header



sub markReadsStatus {
    my ($maxReadsID, $in, $out) = @_;

    my %pairedreads_status = (); # holds seq.xxxx
    my %forwardreads_status = (); # holds seq.xxxa
    my %reversereads_status = (); # holds seq.xxxb


    open IN, $in or croak "Error open $in.\n";
    open OUT, '>'.$out or croak "Error open $out.\n";

    my $counter = 0;

## retrieve read sequence ID
## we assume the read ID is like 12345a, 12345b

while (my $line = <IN>){
    # col 10 is the read ID
    my @a = split/\t/,$line;
    my $id = $a[9];
    $id=~/(\d+)(.)/;  
    print $id,"\t",$1,"\t",$2,"\n";
    $forwardreads_status{$1}++  if $2 eq'a';
    $reversereads_status{$1}++  if $2 eq'b';
    $pairedreads_status{$1}++;

    }

    <STDIN>;
$counter=0;

foreach my $key (1..$maxReadsID){
    if ( !$pairedreads_status{$key} ) {
        $pairedreads_status{$key}='nomap';
    }
    print OUT $key, "\t",$pairedreads_status{$key},"\t";

    if ( $pairedreads_status{$key} == 1 ) {
        if ( $forwardreads_status{$key} == 1 ) {
            $pairedreads_status{$key}="u_forward";
        }else {
            $pairedreads_status{$key}="u_reverse";
        }

    } else {

        if ( $forwardreads_status{$key}==1 && $reversereads_status{$key}==1){
            $pairedreads_status{$key}="u_pair";
        } elsif ($forwardreads_status{$key}>=1 && $reversereads_status{$key}>=1){
            $pairedreads_status{$key}="nu_pair";
        } elsif ($forwardreads_status{$key}>1 && $reversereads_status{$key}<1) {
            $pairedreads_status{$key}="nu_forward";
        } elsif ( $forwardreads_status{$key}<1 && $reversereads_status{$key}>1) {
            $pairedreads_status{$key}="nu_reverse";
        }

    }
    print OUT $pairedreads_status{$key},"\n";


}

    close IN;
    close OUT;
}


### Find status, which is one of 
###  u_pair, u_forward, u_reverse, nu_pair, nu_forward, nu_reverse, nomap
###

sub findStatus {
    my ($status, $in, $out) = @_;

    open IN, $in or die "Error open $in. ";
    open OUT, '>'.$out or die "Error open $out. ";

    while(<IN>){
	my @a=split/\t/;
	if ( $a[2] =~/^$status/){
	    print OUT $_;
	}

    }

    close IN;
    close OUT;
}
