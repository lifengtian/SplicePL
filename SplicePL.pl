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

#use Bio::SeqIO;


require Getopt::Long;
my $USAGE   = '';
my $VERSION = '';
my $FORWARD_FILENAME = undef;
my $REVERSE_FILENAME = undef;
my $GENOME_DIR = undef;
my $GENOME_2BIT_FILENAME = undef; 
my $PROCESSES = 1;
my $TILESIZE = 11;
my $STEPSIZE = 5;
my $MINSCORE = 30;
my $FLANKSIZE = 10;
my $INTRON_MIN = 50;
my $INTRON_MAX = 750000;

my $getopt  = Getopt::Long::GetOptions(	
        'help|usage'    => \$USAGE,
	'forward=s'     => \$FORWARD_FILENAME,
	'reverse=s'     => \$REVERSE_FILENAME,
	'genome=s'      => \$GENOME_DIR,
	'genome_2bit=s' => \$GENOME_2BIT_FILENAME,
	'processes=i'   => \$PROCESSES,
	'tilesize=i'    => \$TILESIZE,
	'stepsize=i'    => \$STEPSIZE,
	'minscore=i'    => \$MINSCORE,
	'flanksize=i'   => \$FLANKSIZE,
	'intron_min=i'  => \$INTRON_MIN,
	'intron_max=i'  => \$INTRON_MAX,
	);




# commandline usage
&usage(), exit(1) if ( $USAGE ) ;


if ( not $GENOME_DIR or not $GENOME_2BIT_FILENAME or not $FORWARD_FILENAME or not $REVERSE_FILENAME ) {
if ( not $GENOME_DIR ) {
    print "Set --genome=/path/to/fasta .\n";
} 

if (  not $GENOME_2BIT_FILENAME ) {
         print "Set --genome_2bit=2bit_file_name \n";
         } 
         
         if (  not $FORWARD_FILENAME ) {
             print "Set --forward=/forward_read_file_name \n";
             } 
             
             if ( not $REVERSE_FILENAME ) {
                 print "--reverse=/reverse_read_file_name \n";
             }
&usage();
exit(1);
}

print "FORWARD_FILENAME = $FORWARD_FILENAME 
REVERSE_FILENAME = $REVERSE_FILENAME
GENOME_DIR = $GENOME_DIR 
GENOME_2BIT_FILENAME = $GENOME_2BIT_FILENAME 
PROCESSES = $PROCESSES
TILESIZE = $TILESIZE
STEPSIZE = $STEPSIZE
MINSCORE = $MINSCORE
FLANKSIZE = $FLANKSIZE
INTRON_MIN = $INTRON_MIN
INTRON_MAX = $INTRON_MAX
";

### let's check genome_2bit size and total available physical memory
if ( -f $GENOME_DIR.'/'.$GENOME_2BIT_FILENAME ) {
    my $filesize = -s $GENOME_DIR.'/'.$GENOME_2BIT_FILENAME ;
    check_total_memory();
} else {
        croak "$GENOME_DIR.'/'.$GENOME_2BIT_FILENAME does NOT exist. Exit!\n";
}

exit();
############################### Check Executables, dababase etc. #############################

my $current_path = abs_path();
my $temp         = $current_path . '/temp';
my $fasta_fn    = $temp.'/pairedreads.fa';
my $psl_fn       = $temp . '/pairedreads.psl';

if ( !-d $temp ) {
    mkdir($temp);
}

exit(0);

### have to fix is_fasta and prepare_pairedreads next!!!
if ( is_fasta($current_path.'/'.$FORWARD_FILENAME,  $current_path.'/'.$REVERSE_FILENAME) ) {
       &prepare_pairedreads($current_path.'/'.$FORWARD_FILENAME,  $current_path.'/'.$REVERSE_FILENAME, $fasta_fn);
} else {
    croak "Not fasta file!\n";
}



&split_fasta( $fasta_fn, $PROCESSES, $current_path );

&run_blat( $fasta_fn, $GENOME_DIR.'/'.$GENOME_2BIT_FILENAME, $PROCESSES, $TILESIZE, $STEPSIZE, $MINSCORE,
    $current_path );

&mark_reads_status( $psl_fn, $fasta_fn );

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

my $junction_path = $current_path . '/junctions';
mkdir($junction_path);
&find_junctions( $temp . '/u_pair_ALL.psl',
    $junction_path . '/junctions.junc' , $INTRON_MIN, $INTRON_MAX, $FLANKSIZE );
    
&unique_junctions( $junction_path . '/junctions.junc',
    $junction_path . '/unique.junc' );

&filter_splicesite( $junction_path . '/unique.junc',
    $junction_path . '/unique_GTAG.junc' );

#################################### Subroutines  #################################

#  = head2 Title : Usage : Function : Returns : Args :=cut

=head2   run_blat
 Title   : run_blat
 Usage   : run_blat("/path/to/run_folder",6)
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

## You must run SplicePL in the same folder as fasta file

    unless ( -f $cp . '/' . $fasta ) {
        croak "Can't find $fasta file in $cp. Exit!";
    }

    if ( !-d $cp . '/scripts' ) {
        mkdir( $cp . '/scripts' );
    }

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

sub split_fasta {
    my ( $fasta, $chunk, $cp ) = @_;

    open IN, $cp . '/' . $fasta or croak "Error open $fasta \n";
    my $fasta_no = 0;
    map { $fasta_no++ if /^>/ } <IN>;
    close IN;
    print "No of fasta records: ", $fasta_no, "\n";

    open IN, $cp . '/' . $fasta or croak "Error open $fasta \n";
    my $chunk_size = int( $fasta_no / $chunk );

    croak "Can't split $fasta_no to less than 2 chunks. " if $chunk <= 1;
    $chunk_size++ if $chunk_size % 2 == 1;

    foreach my $i ( 1 .. ( $chunk - 1 ) ) {
        my $out = $cp . '/' . $fasta . "." . $i;
        open OUT, '>' . $out or croak "Error open $out \n";

        foreach my $j ( 1 .. $chunk_size ) {
            my $l = <IN>;
            print OUT $l;
            $l = <IN>;
            print OUT $l;
        }
        close(OUT);
    }

    my $out = $cp . '/' . $fasta . "." . $chunk;
    open OUT, '>' . $out or croak "Error open $out \n";

    map { print OUT } <IN>;

    close(OUT);
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

=head2 
Title :  
Usage : 
Function : find junctions from a psl file 
 col  9: strand
 col 10: seqID
 col 11: qSize
 col 14: chr
 col 18: block no
 col 19: blockSize
 col 21: tStart  (target start position)
 col start from 1 (not zero)
 request blocksize at least 10 bases and gap between blocks at least 20 bases for long reads
 may change gap size to 50bp in the future
Returns : 
Args :
=cut

sub find_junctions {
    my ( $psl, $out, $mingap, $maxgap, $minblocksize ) = @_;

    my $first_row = 1;

    open IN,  $psl       or croak "Error open $psl \n";
    open OUT, '>' . $out or croak "Error open $out \n";

    
    while (<IN>) {
        if ( $first_row == 1 ) {
            my @a = split /\t/;
            if ( $a[10] <= 60 && $a[10] >= 20 ) {
                $minblocksize = 10 if $minblocksize < 20;
            }
            elsif ( $a[10] > 60 )
            {    #query size >=60, might need larger blocksize
                $minblocksize = 20 if $minblocksize < 20;
            }
            else {
                die "PSL file has abnormal qSize: $a[10]. Quit!";
            }
            print "minblocksize=$minblocksize\n";
            $first_row = undef;
        }
        my @a = split /\t/;
        my @tstarts = split /,/, $a[20];

        my @blocksize = split /,/, $a[18];
        my $blockno = $a[17];
        if ( $blockno > 1 ) {
            foreach ( 0 .. $#tstarts - 1 ) {
                my $gapsize =
                  $tstarts[ $_ + 1 ] + 1 - ( $tstarts[$_] + $blocksize[$_] );
                if (   $gapsize >= $mingap
                    && $gapsize <= $maxgap
                    && $blocksize[$_] >= $minblocksize
                    && $blocksize[ $_ + 1 ] >= $minblocksize )
                {

                    print OUT $a[13], ":", $tstarts[$_] + $blocksize[$_], "-",
                      $tstarts[ $_ + 1 ] + 1, "\t", $a[9], "\t", $a[8], "\t",
                      $tstarts[ $_ + 1 ] + 1 -
                      ( $tstarts[$_] + $blocksize[$_] ),
                      "\n";
                }
            }
        }
    }

}

=head2 unique_junctions 
Title : unique_junctions
Usage : 
Function : 
Returns : 
Args :
=cut

sub unique_junctions {
    my ( $in, $out ) = @_;

    my %h;
    my %line;
    my %l;

    open IN,  $in        or croak "Error open $in \n";
    open OUT, '>' . $out or croak "Error open $out \n";

    while (<IN>) {
        if (/^chr/) {
            chomp;
            my @a = split /\t/;
            $h{ $a[0] }++;
            $line{ $a[0] } = $_;
            $l{ $a[0] } .= '_' . $a[1];
        }
    }

    foreach ( sort keys %h ) {
        if ( $h{$_} ) {
            print OUT $line{$_}, "\t", $l{$_}, "\n";
        }
    }

    close IN;
    close OUT;
}

### Filter the junctions with the splicing site patterns
### 1. GT...AG
### 2. GC...AG
### 3. AT...AC
###

sub filter_splicesite {
    my ( $junc, $out, $dir ) = @_;
### Global path

    my ( $splus, $sneg );

    open IN,  $junc      or croak "Error open $junc \n";
    open OUT, '>' . $out or croak "Error open $out \n";


    my @names = (
        1,  2,  3,  4,  5,  6,  7,  8,  9,   10,  11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y', 'M'
    );

    my %seq_obj = undef ;    # store the reference to sequence objects

### Create Seq object for each chromosome
    foreach my $chr (@names) {
        my $file = $dir . '/chr' . $chr . '.fa';

        print "Loading chromosome $file\n";
        open CHR, $file or croak "Error open $file \n";
        my $id = <CHR>;
        croak "First line of fasta file $file doesn't start with >\n" unless ( $id =~/^>/ );
        my $seq = '';
        while ( <CHR> ){
            chomp;
            $seq .= $_;
        }
        $seq_obj{ 'chr' . $chr } = $seq ;
    }

### We were not analyzing strand-specific RNA-seq data, so
### have to search both plus and minus strand for
### splicing site patterns

    while (<IN>) {
        chomp;
        my @a = split /\s+/;
        my ( $chr, $start, $end ) = split /[:-]/, $a[0];
        my $strand = $a[2];
        my $s;
        $splus =
            substr( $seq_obj{$chr}, $start, 2 )
          . substr( $seq_obj{$chr}, $end - 3, 2 );
        $sneg = &revcom($splus);

        #print $splus,"\n",$sneg;<STDIN>;
        #print STDERR ">",$a[0],"\n",$splus,"\n";

        if ( $splus =~ /^GTAG$/i ) {
            print OUT $_, "\tGOODSITEPLUS" . $strand . "GTAGSITE\n";
        }
        elsif ( $sneg =~ /^GTAG$/i ) {
            print OUT $_, "\tGOODSITENEG" . $strand . "GTAGSITE\n";
        }
        elsif ( $splus =~ /^GCAG$/i ) {
            print OUT $_, "\tGOODSITEPLUS" . $strand . "GCAGSITE\n";
        }
        elsif ( $sneg =~ /^GCAG$/i ) {
            print OUT $_, "\tGOODSITENEG" . $strand . "GCAGSITE\n";
        }
        elsif ( $splus =~ /^ATAC$/i ) {
            print OUT $_, "\tGOODSITEPLUS" . $strand . "ATACSITE\n";
        }
        elsif ( $sneg =~ /^ATAC$/i ) {
            print OUT $_, "\tGOODSITENEG" . $strand . "ATACSITE\n";
        }
        else {
            print OUT $_, "\tBADSITE$strand\n";
        }
    }

    close IN;
    close OUT;
}

sub revcom {
    my ($seq) = @_;

    my $r = reverse($seq);
    $r =~ tr/atcgATCG/tagcTAGC/;

    return $r;

}


sub prepare_pairedreads {
        my ( $f, $r, $out ) = @_;
        
        open F, $f or croak "Error open $f $!";
    
    #check fasta
    #combine each line of fasta to a string, get rid of newline, check characters etc.
    #renumber the fasta record ID, make it like 232f followed by 232r
    #output it to a file pairedreads.fa
    #DONE
    
}


# Retrieves the memory installed on this machine
sub check_total_memory {
    
sub memerror {
		croak "Unable to determine total memory\n";
	}
	

    my ($physical_memory,$swap_memory,$duflags);
	my $os = `uname`;

		if ($os =~ /Linux/) {
			$physical_memory = `free -b | grep Mem | awk '{print \$2}'` or memerror;
		} elsif ($os =~ /Darwin/) {
			$physical_memory = `sysctl -n hw.memsize` or memerror;
		} elsif ($os =~ /NetBSD|OpenBSD/) {
			$physical_memory = `sysctl -n hw.physmem` or memerror;
			if ($physical_memory < 0) {
				$physical_memory = `sysctl -n hw.physmem64` or memerror;
			}
			
		} elsif ($os =~ /BSD/) {
			$physical_memory = `sysctl -n hw.realmem`;
			
		} elsif ($os =~ /SunOS/) {
			$physical_memory = `/usr/sbin/prtconf | grep Memory | cut -f 3 -d ' '` or memerror;
			chomp($physical_memory);
			$physical_memory = $physical_memory*1024*1024;
		}
	
	chomp($physical_memory);
	print "Total physical memory is ", int ($physical_memory/1024/1024/1024 ), "GB\n";
}


	
sub usage {
	print <<"END_USAGE";
Usage: $0 

--help                  Shows this help message
--forward=name          Name of forward read file, fasta file only
--reverse=name          Name of reverse read file, fasta file only
--genome=dir            Directory of genome sequences for individual chromosome, ie, chr1.fa, chr2.fa, ..., chr22.fa, chrX.fa, chrY.fa, chrM.fa
--genome_2bit=name      Name of genome sequence file in 2bit format
--processes             Number of processes for BLAT (default is 1)
--tilesize              Number of tileSize (default is 11)
--stepsize              Number of stepSize (default is 5)
--minscore              Number of minScore (default is 30)
--flanksize             Minimum length of splice junction flanking sequence (default is 10)
--intron_min            Minimum length of intron (default is 50)
--intron_max            Maximum length of intron (default is 750000)

END_USAGE
}