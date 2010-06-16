#!/usr/bin/perl 

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

package SplicePL;

use strict;
use Carp;
use Cwd 'abs_path';
use Env '@PATH';

my $usage =
"perl /path/to/SplicePL.pl <fasta file name> </path/to/database> <no. of threads> <tileSize> <stepSize> <minScore>\n";
if ( @ARGV != 6 ) {
    die $usage;
}

my $fasta_fn = $ARGV[0];    # input fasta file
my $db       = $ARGV[1];
my $chunk    = $ARGV[2];    # number of threads to parallelize
my $tilesize = $ARGV[3];    # BLAT tileSize
my $stepsize = $ARGV[4];    # BLAT stepSize
my $minscore = $ARGV[5];    # BLAT minScore

############################### Check Executables, dababase etc. #############################

my $current_path = abs_path();
my $temp         = $current_path . '/temp';
my $psl_fn       = $temp . '/' . $fasta_fn . '.psl';

if ( !-d $temp ) {
    mkdir($temp);
}

=head1
&breakup_fasta( $fasta_fn, $chunk, $current_path );

&run_blat( $fasta_fn, $db, $chunk, $tilesize, $stepsize, $minscore,
    $current_path );

&mark_reads_status( $psl_fn, $fasta_fn );

=cut

## unique-reads
&find_status( 'u_pair', $psl_fn, $temp . '/u_pair' );
&retrieve_status_from_psl( $temp . '/u_pair', $psl_fn, $temp . '/u_pair.psl' );
&filter_paired_reads(
    $temp . '/u_pair.psl',
    $temp . '/u_pair_consistent.psl',
    $temp . '/u_pair_inconsistent'
);

## non-unique reads
&find_status( 'nu_pair', $psl_fn, $temp . '/nu_pair' );
&retrieve_status_from_psl( $temp . '/nu_pair', $psl_fn,
    $temp . '/nu_pair.psl' );
&filter_nu_reads_from_psl(
    $temp . '/nu_pair.psl',
    $temp . '/nu_pair_tophit.psl',
    $temp . '/nu_pair_tophit_excluded.psl'
);
### find unique reads after filtering nu_reads
&mark_reads_status( $temp . '/nu_pair_tophit.psl', $fasta_fn );

&find_status(
    'u_pair',
    $temp . '/nu_pair_tophit.psl',
    $temp . '/nu_pair_tophit_u_pair'
);

&retrieve_status_from_psl(
    $temp . '/nu_pair_tophit_u_pair',
    $temp . '/nu_pair_tophit.psl',
    $temp . '/u_pair_from_nu_pair_tophit.psl'
);

&filter_paired_reads(
    $temp . '/u_pair_from_nu_pair_tophit.psl',
    $temp . '/u_pair_from_nu_pair_tophit_consistent.psl',
    $temp . '/u_pair_from_nu_pair_tophit_inconsistent'
);

system( 'cat ' 
      . $temp
      . '/u_pair_consistent.psl '
      . $temp
      . '/u_pair_from_nu_pair_tophit_consistent.psl >'
      . $temp
      . '/u_pair_ALL.psl' );

## find junctions

my $junction_path = $current_path.'/junctions';
mkdir ( $junction_path );
&find_junctions( $temp . '/u_pair_ALL.psl', $junction_path.'/junctions.junc'  );
&unique_junctions( $junction_path.'/junctions.junc', $junction_path.'/unique.junc' );
#&filter_splicesite ($junction_path.'/unique.junc');

=head1

  $jp / junctions $uniquejunctions $jp / junctions >
  $jp / uniquejunctions $filter_GTAG2 $jp / uniquejunctions |
  grep -h GOOD | awk '{if ($4>50) print $0}' >
  $jp / filtered_junctions date echo Done
=cut


#################################### Subroutines  #################################

#  = head2 Title : Usage : Function : Returns : Args :=cut

=head2   runBlat
 Title   : runBlat
 Usage   : runBlat("/path/to/run_folder",6)
 Function: run multiple Blat processes
 Returns : 
 Args    : First is the directory
           Second is the number of subprocesses
=cut

  sub run_blat {

    my ( $fasta, $db, $chunk, $tilesize, $stepsize, $minscore, $cp ) = @_;

    my $blat_bin = 'blat';

    unless ( grep -x "$_/$blat_bin", @PATH ) {
        croak "Can't find blat executables in @PATH\n";
    }

######  check the databases for BLAT
    unless ( -f $db ) {
        croak "Can't find BLAT database file as $db\n";
    }

######  current path #######
## You must run SplicePL in the same folder as fasta file

    unless ( -f $cp . '/' . $fasta ) {
        croak "Can't find $fasta file in $cp. Exit!";
    }

    if ( !-d $cp . '/scripts' ) {
        mkdir( $cp . '/scripts' );
    }

################################# BLAT parameters #############################################
    my $repmatch = int( 1024.0 * $tilesize / $stepsize );
    my $parameters =
"-noHead -repMatch=$repmatch -minScore=$minscore -stepSize=$stepsize -tileSize=$tilesize";

    for ( my $i = 1 ; $i <= $chunk ; $i++ ) {
        my $command = <<EOF;
#!/bin/sh' 
date >> $cp/run.$i.log
echo blat $db $parameters $cp/$fasta.$i $cp/$fasta.$i.out >> $cp/run.$i.log
blat $db $parameters $cp/$fasta.$i $cp/$fasta.$i.out
date >> $cp/run.$i.log
EOF

        my $of = 'run.' . $i . '.sh';
        open OUT, ">$cp/scripts/" . $of or die "Error open $of for writing. ";
        print OUT $command;
    }

    for my $p ( 1 .. $chunk ) {
        my $pid = fork();
        if ( $pid == -1 ) {
            die;
        }
        elsif ( $pid == 0 ) {

            my $blat_process = "bash $cp/scripts/run.$p.sh ";
            exec $blat_process or die;
        }
    }
    while ( wait() != -1 ) { }
    print "All $chunk BLAT jobs Done \n";

## merge BLAT output files into one file
    my $s    = 'cat ';
    my $temp = $cp . '/temp';

    for my $j ( 1 .. $chunk ) {
        $s .= $cp . '/' . $fasta . '.' . $j . '.out' . ' ';
    }

    my $psl_fn   = $temp . '/' . $fasta . '.psl';
    my $comm_cat = $s . '> ' . $psl_fn;
    print $comm_cat, "\n";
    system($comm_cat);

}

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
    my ( $fastafile, $numpieces, $cp ) = @_;
    open( INFILE, $cp . '/' . $fastafile );

    print STDERR "counting fasta records in the $cp/$fastafile\n";

    my $filesize = `grep '>' $cp/$fastafile | wc -l`;
    chomp($filesize);
    $filesize =~ s/\s*//;
    print STDERR "No of fasta records: $filesize\n";

    my $numseqs   = $filesize;
    my $piecesize = int( $numseqs / $numpieces );

    if ( $piecesize < 1 ) {
        die "Quit! Not enough pieces to break into $numpieces!";
    }

    print STDERR
      "processing in $numpieces pieces of approx $piecesize size each.\n";
    if ( $piecesize % 2 == 1 ) {
        $piecesize++;
    }
    my $bflag = 0;
    for ( my $i = 1 ; $i < $numpieces ; $i++ ) {
        my $outfilename = $cp . '/' . $fastafile . "." . $i;
        open( OUTFILE, ">$outfilename" );
        for ( my $j = 0 ; $j < $piecesize ; $j++ ) {
            my $line = <INFILE>;
            print OUTFILE $line;
            $line = <INFILE>;
            print OUTFILE $line;
        }
        close(OUTFILE);
    }
    my $outfilename = $cp . '/' . $fastafile . "." . $numpieces;
    open( OUTFILE, ">$outfilename" );
    while ( my $line2 = <INFILE> ) {
        print OUTFILE $line2;
    }
    close(OUTFILE);
    return 0;
}

=head2 mark_reads_status
 Title   : mark_reads_status
 Usage   : mark_reads_status("/path/to/psl_file","/path/to/fasta_file")
 Function: 
           Each read will be marked as
            u_pair
            u_forward
            u_reverse
            nu_pair
            nu_forward
            nu_reverse
            nomap
 Returns : 
 Args    :  First string is the path to input PSL file 
            Second is the path to the fasta file
=cut

sub mark_reads_status {
    my ( $psl, $fasta ) = @_;

    my $status = $psl . '.status';

    my %pairedreads_status  = ();
    my %forwardreads_status = ();
    my %reversereads_status = ();

    open PSL,   $psl          or croak "Error open $psl.\n";
    open OUT,   '>' . $status or croak "Error open $status.\n";
    open FASTA, $fasta        or croak "Error open $fasta.\n";

    my $counter       = 0;
    my $maxFastaID    = 0;    # in fasta file
    my $fastaTotalSeq = 0;    # total number of sequences in fasta file
    my $maxReadID     = 0;    # in psl file

    print join( "\t", $psl, $status, $fasta ), "\n";

## retrieve read sequence ID
## we assume the read ID is like 12345a, 12345b and it is CONTINUALLY numbered from 1 to maxReadID
    while (<FASTA>) {
        if (/^>\D*(\d+)/) {
            $fastaTotalSeq++;
            if ( $maxFastaID < $1 ) {
                $maxFastaID = $1;
            }
        }
    }

    print "total fasta sequences: ", $fastaTotalSeq, " maxFastaID: ",
      $maxFastaID, "\n";

    while (<PSL>) {

        # col 10 is the read ID
        my $id = ( split /\t/ )[9];

        if ( $id =~ /(\d+)(.)/ ) {

            #print $id, "\t", $1, "\t", $2, "\n";
            if ( $1 > $maxReadID ) {
                $maxReadID = $1;
            }
            if ( $2 eq 'f' ) {
                $forwardreads_status{$1}++;
            }
            elsif ( $2 eq 'r' ) {
                $reversereads_status{$1}++;
            }
            $pairedreads_status{$1}++;
        }
    }

    #print "maxReadID: $maxReadID \n";

    if ( $maxReadID > $maxFastaID ) {
        croak "maxReadID larger than maxFastaID. Exit!\n";
    }

    foreach my $key ( 1 .. $maxFastaID ) {
        if ( !$pairedreads_status{$key} ) {
            $pairedreads_status{$key} = "nomap";
        }
        print OUT $key, "\t", $pairedreads_status{$key}, "\t";

        if ( $pairedreads_status{$key} == 1 ) {
            if ( $forwardreads_status{$key} == 1 ) {
                $pairedreads_status{$key} = "u_forward";
            }
            else {
                $pairedreads_status{$key} = "u_reverse";
            }

        }
        else {

            if (   $forwardreads_status{$key} == 1
                && $reversereads_status{$key} == 1 )
            {
                $pairedreads_status{$key} = "u_pair";
            }
            elsif ($forwardreads_status{$key} >= 1
                && $reversereads_status{$key} >= 1 )
            {
                $pairedreads_status{$key} = "nu_pair";
            }
            elsif ($forwardreads_status{$key} > 1
                && $reversereads_status{$key} < 1 )
            {
                $pairedreads_status{$key} = "nu_forward";
            }
            elsif ($forwardreads_status{$key} < 1
                && $reversereads_status{$key} > 1 )
            {
                $pairedreads_status{$key} = "nu_reverse";
            }

        }

        print OUT $pairedreads_status{$key}, "\n";

    }

    close IN;
    close OUT;
}

### Find status, which is one of
###  u_pair, u_forward, u_reverse, nu_pair, nu_forward, nu_reverse, nomap
###

sub find_status {
    my ( $status, $psl, $out ) = @_;

    my $status_psl = $psl . '.status';
    open IN,  $status_psl or croak "Error open $status_psl. ";
    open OUT, '>' . $out  or croak "Error open $out. ";

    while (<IN>) {
        my $a = ( split /\t/ )[2];
        if ( $a =~ /^$status/ ) {
            print OUT $_;
        }

    }

    close IN;
    close OUT;
}

sub retrieve_status_from_psl {
### given a list of u_pair read ID such as 123
### get the blat result in the psl file
### psl file has no HEADERS by -noHead when running blat

    my ( $in, $psl, $out ) = @_;
    my %ids = ();

    open PSL, $psl       or croak "Error open $psl \n";
    open ID,  $in        or croak "Error open $in \n";
    open OUT, '>' . $out or croak "Error open $out \n";

    while (<ID>) {
        chomp;
        my $id = ( split /\t/ )[0];

        my $idf = $id . 'f';
        my $idr = $id . 'r';
        $ids{$idf}++;
        $ids{$idr}++;
    }

    while (<PSL>) {
        my $i = ( split /\t/ )[9];
        if ( $ids{$i} ) {
            print OUT;
        }
    }

    close ID;
    close PSL;
    close OUT;

}

## exclude reads
##   1. mapped to different chromosome
##   2. mapped to same strand
##   3. too far away

sub filter_paired_reads {
    my ( $psl, $out1, $out2 ) = @_;

    my $maxdist = 750000;

    print "filter_paired_reads $psl, $out1, $out2 \n";

    open IN,   $psl        or croak "Error open $psl \n";
    open OUT1, '>' . $out1 or croak "Error open $out1 \n";
    open OUT2, '>' . $out2 or croak "Error open $out2 \n";

    my $count = 1;

    while ( my $line1 = <IN> ) {
        my ( $for_id, $for_str, $for_chr, $for_start ) =
          ( split /\t/, $line1 )[ 9, 8, 13, 15 ];

        my $m1;
        my $m2;

        my $line2 = <IN>;
        my ( $rev_id, $rev_str, $rev_chr, $rev_start ) =
          ( split /\t/, $line2 )[ 9, 8, 13, 15 ];

        if ( $for_id =~ /(\d+)f/ ) {
            $m1 = $1;
        }
        else {
            croak $for_id,
              "does not match the format of forward ID (e.g., 1234f) \n";
        }

        if ( $rev_id =~ /(\d+)r/ ) {
            $m2 = $1;
        }
        else {
            croak $rev_id,
              "does not match the format of reverse ID (e.g., 1234r) \n";
        }

        if ( $m1 ne $m2 ) {
            croak "forward read ID does not match reverse read ID. Exit.\n";
        }

 #    print join("\t",$for_id,$for_str,$for_chr,$rev_id,$rev_str,$rev_chr),"\n";

        if ( $for_str eq $rev_str || $for_chr ne $rev_chr ) {
            ## strand and chromosome
            if ( $for_str eq $rev_str && $for_chr ne $rev_chr ) {
                print OUT2 $count, " same str: ",
                  join( "\t", $for_id, $for_str, $rev_id, $rev_str ), "\t";
                print OUT2 $count++, " diff chr: ",
                  join( "\t", $for_id, $for_chr, $rev_id, $rev_chr ), "\n";
            }
            elsif ( $for_str eq $rev_str ) {
                print OUT2 $count++, " same str: ",
                  join( "\t", $for_id, $for_str, $rev_id, $rev_str ), "\n";
            }
            elsif ( $for_chr ne $rev_chr ) {
                print OUT2 $count++, " diff chr: ",
                  join( "\t", $for_id, $for_chr, $rev_id, $rev_chr ), "\n";
            }

        }
        else {
            ## distance
            if ( abs( $for_start - $rev_start ) >= $maxdist ) {
                print OUT2 $count++, " too far: ",
                  join( "\t", $for_id, $for_chr, $rev_id, $rev_chr ), "\n";
            }
            else {
                print OUT1 $line1, $line2;
            }
        }
    }

}

sub filter_nu_reads_from_psl {
    my ( $psl, $out1, $out2 ) = @_;
    my %h;
    my %count_max;
    my $previous;
    my $maxline;

    open IN,   $psl        or croak "Error open $psl \n";
    open OUT1, '>' . $out1 or croak "Error open $out1 \n";
    open OUT2, '>' . $out2 or croak "Error open $out2 \n";

    my @a = split /\t/, <IN>;
    $h{ $a[9] } = $a[0];    ## col 10 is the seqID  ## col1 is the score
    $previous = $a[9];

    while (<IN>) {
        my @a = split /\t/;
        if ( !$h{ $a[9] } ) {
            if ( $count_max{$previous} ) {
                print OUT2 '?' . $previous, "\t", $count_max{$previous};
            }
            else {
                print OUT1 $maxline;
            }
            $previous = $a[9];
        }
        if ( $h{ $a[9] } < $a[0] ) {
            $h{ $a[9] }         = $a[0];
            $maxline            = $_;
            $count_max{ $a[9] } = undef;
        }
        elsif ( $h{ $a[9] } == $a[0] ) {
            $count_max{ $a[9] } = $maxline . $a[9] . "\t" . $_;
        }
    }

}



## find junctions from a psl file
## col  9: strand
## col 10: seqID
## col 11: qSize
## col 14: chr
## col 18: block no
## col 19: blockSize
## col 21: tStart  (target start position)
## col start from 1 (not zero)
## request blocksize at least 10 bases and gap between blocks at least 20 bases for long reads
## may change gap size to 50bp in the future


sub find_junctions {
    my ($psl, $out) = @_;


## constant ###
    my $minGap       = 50;       ### in the final analysis, a gap of 50 is used.
    my $maxGap       = 750000;
    my $minBlockSize = 10;

    my $first_row = 1;

    open IN, $psl or croak "Error open $psl \n";
    open OUT, '>'.$out or croak "Error open $out \n";

    while (<IN>) {
        if ( $first_row == 1 ) {
            my @a = split /\t/;
            if ( $a[10] <= 60 && $a[10] >= 20 ) {
                $minBlockSize = 10;
            }
            elsif ( $a[10] > 60 ) {
                $minBlockSize = 10;
            }
            else {
                die "PSL file has abnormal qSize: $a[10]. Quit!";
            }
            print STDERR "minBlockSize=$minBlockSize\n";
            $first_row = undef;
        }
        my @a = split /\t/;
        my @tstarts = split /,/, $a[20];

        my @blockSize = split /,/, $a[18];
        my $blockno = $a[17];
        if ( $blockno > 1 ) {
            foreach ( 0 .. $#tstarts - 1 ) {
                my $gapsize =
                  $tstarts[ $_ + 1 ] + 1 - ( $tstarts[$_] + $blockSize[$_] );
                print OUT $a[13], ":", $tstarts[$_] + $blockSize[$_], "-",
                  $tstarts[ $_ + 1 ] + 1, "\t", $a[9], "\t", $a[8], "\t",
                  $tstarts[ $_ + 1 ] + 1 - ( $tstarts[$_] + $blockSize[$_] ),
                  "\n"
                  if ( $gapsize >= $minGap
                    && $gapsize <= $maxGap
                    && $blockSize[$_] >= $minBlockSize
                    && $blockSize[ $_ + 1 ] >= $minBlockSize );
            }
        }
    }

}


sub unique_junctions {
    my ($in, $out ) = @_;
    
    my %h;
    my %line;
    my %l;

    open IN, $in or croak "Error open $in \n";
    open OUT, '>'.$out or croak "Error open $out \n";

while(<IN>){
if (/^chr/ ){
        chomp;
        my @a=split/\t/;
        $h{$a[0]}++;
        $line{$a[0]}=$_;
        $l{$a[0]}.='_'.$a[1];
}
}


foreach (sort keys %h){
        if ($h{$_} ) {
                print OUT $line{$_},"\t",$l{$_},"\n";
        }
}

    close IN;
    close OUT;
}
