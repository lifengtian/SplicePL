#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# PERL package for SplicePL
#
# Please direct questions and support issues to <lifeng2@mail.med.upenn.edu>
#
# Copyright Lifeng Tian
#
# You may distribute this module under the same terms as PERL itself

# POD documentation - main docs before the code

=head1 NAME

=head1 DESCRIPTION


=head1 Check List 
Before running the SplicePL, make sure:
1. BLAT binary is in the PATH environment variable.
2. SplicePL is in the PATH, and executable!
3. BLAT database file, i.e., the xxx.2bit file

Users need to make sure that several binaries are installed, including BLAT, 
databases including
1. genome/chr1 to chrM.fa
2. genome/genome.2bit

=cut 

package SplicePL;

use strict;
use Carp;
use Cwd 'abs_path';
use Env '@PATH';
use FileHandle;
use POSIX qw(strftime);

use Getopt::Long;

### run() was called only when user execute "perl SplicePL.pm "
&run() unless caller();

sub run
{
    my $USAGE                = '';
    my $VERSION              = '1.0';
    my $FORWARD_FILENAME     = undef;
    my $REVERSE_FILENAME     = undef;
    my $GENOME_DIR           = undef;
    my $GENOME_2BIT_FILENAME = undef;
    my $PROCESSES            = undef;
    my $TILESIZE             = 11;
    my $STEPSIZE             = 5;
    my $MINSCORE             = 30;
    my $FLANKSIZE            = 10;
    my $INTRON_MIN           = 50;
    my $INTRON_MAX           = 750000;
    my $REPMATCH             = undef;

    my $getopt =
      Getopt::Long::GetOptions(
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
                               'repmatch'      => \$REPMATCH
                              );

    # commandline usage
    &usage(), exit(1) if ($USAGE);

    #
    if (   not $GENOME_DIR
        or not $GENOME_2BIT_FILENAME
        or not $FORWARD_FILENAME
        or not $REVERSE_FILENAME)
    {
        if (not $GENOME_DIR)
        {
            print "Set --genome=/path/to/fasta .\n";
        }

        if (not $GENOME_2BIT_FILENAME)
        {
            print "Set --genome_2bit=2bit_file_name \n";
        }

        if (not $FORWARD_FILENAME)
        {
            print "Set --forward=/forward_read_file_name \n";
        }

        if (not $REVERSE_FILENAME)
        {
            print "--reverse=/reverse_read_file_name \n";
        }
        &usage();
        exit(1);
    }

    my $start_time = time();
    &mlog("SplicePL($VERSION) started");

    print "Parameters\nFORWARD_FILENAME = $FORWARD_FILENAME 
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
\n";

    my $error = undef;

### check individual human genome chromosome fasta file
    my @chr_names = (
                     1,  2,  3,   4,   5,  6,  7,  8,  9,  10,
                     11, 12, 13,  14,  15, 16, 17, 18, 19, 20,
                     21, 22, 'X', 'Y', 'M'
                    );
    foreach my $chr (@chr_names)
    {
        $error .= "!!ERROR!! chr$chr does not exist in $GENOME_DIR.\n"
          unless -f $GENOME_DIR . '/chr' . $chr . '.fa';
    }

    ### check executable for BLAT
    my $blat_bin = 'blat';
    unless (grep -x "$_/$blat_bin", @PATH)
    {
        $error .= "Can't find blat executables in @PATH\n";
    }

### check genome_2bit size and total available physical memory
    if (-f $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME)
    {
        my $filesize = -s $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME;
        mlog("== MEMORY ==", $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME,
             " size ", int($filesize / 1024 / 1024), "MB");
        my $total_mem       = int(&check_total_memory() / 1024 / 1024 / 1024);
        my $one_process_mem = int($filesize * 6 / 1024 / 1024 / 1024);
        mlog("== MEMORY == Total physical memory is $total_mem GB");
        if ($total_mem <= 3)
        {
            mlog(
                "== MEMORY == Total physical memory < 3G , not enough to load the genome. Exit"
            );
            exit(1);
        }
        if (!$PROCESSES)
        {
            $PROCESSES = int($total_mem / $one_process_mem);
            mlog(
                "== MEMORY == $PROCESSES processes will be used in the analysis"
            );
        }
        mlog("== MEMORY == Predicted memory usage is ",
             $PROCESSES * $one_process_mem, "GB");

        if ($total_mem <= $PROCESSES * $one_process_mem)
        {
            mlog(
                "== MEMORY == Total physical memory is not enough to run $PROCESSES processes simultaneously. 
    The maximum number of processes you can run is ",
                int($total_mem / $one_process_mem)
                );
            exit(1);
        }

    }
    else
    {
        $error .=
          "!!ERROR!! $GENOME_DIR.'/'.$GENOME_2BIT_FILENAME does NOT exist.\n";
    }

    if ($error)
    {
        mlog($error);
        exit(1);
    }

### Create temporary folder for intermediate files

    my $current_path = abs_path();
    my $temp         = $current_path . '/temp';
    my $fasta_fn     = $current_path . '/pairedreads.fa';
    my $psl_fn       = $temp . '/pairedreads.psl';

    if (!-d $temp)
    {
        mkdir($temp);
    }
    else
    {
        mlog("== CREATE TEMP FOLDER == $temp exists");
    }

### Start the pipeline here
    &mcall(
           "PREPARE PARED READS",
           \&prepare_paired_reads,
           $current_path . '/' . $FORWARD_FILENAME,
           $current_path . '/' . $REVERSE_FILENAME,
           $fasta_fn
          );

    &mcall("SPLIT FASTA", \&split_fasta, $fasta_fn, $PROCESSES, $current_path);

    &mcall(
           "RUN BLAT", \&run_blat,
           $fasta_fn,  $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME,
           $PROCESSES, $TILESIZE,
           $STEPSIZE,  $MINSCORE,
           $REPMATCH,  $current_path
          );

    &mcall("MARK READS ALIGNMENT STATUS",
           \&mark_reads_status, $psl_fn, $fasta_fn);

## unique-reads
    &mcall("FIND UNIQUE PAIRED READS",
           \&find_status, 'u_pair', $psl_fn, $temp . '/u_pair');
    &mcall("RETRIEVE UNIQUE ALIGNMENTS",
           \&retrieve_status_from_psl, $temp . '/u_pair',
           $psl_fn, $temp . '/u_pair.psl');
    &mcall(
           "FILTER INCONSISTENT PAIRED READS",
           \&filter_paired_reads,
           $temp . '/u_pair.psl',
           $temp . '/u_pair_consistent.psl',
           $temp . '/u_pair_inconsistent',
           $INTRON_MAX
          );

## non-unique reads
    &mcall("FIND NON-UNIQUE PAIRED READS",
           \&find_status, 'nu_pair', $psl_fn, $temp . '/nu_pair');
    &mcall("RETRIEVE NON-UNIQUE ALIGNMENTS",
           \&retrieve_status_from_psl, $temp . '/nu_pair',
           $psl_fn, $temp . '/nu_pair.psl');
    &mcall(
           "FIND THE BEST LOCATION FOR NON-UNIQUE READS",
           \&filter_nu_reads_from_psl,
           $temp . '/nu_pair.psl',
           $temp . '/nu_pair_tophit.psl',
           $temp . '/nu_pair_tophit_excluded.psl'
          );
### find unique reads after filtering nu_reads
    &mcall("MARK READS ALIGNMENT STATUS", \&mark_reads_status,
           $temp . '/nu_pair_tophit.psl', $fasta_fn);

    &mcall("FIND UNIQUE PAIRED READS",
           \&find_status, 'u_pair',
           $temp . '/nu_pair_tophit.psl',
           $temp . '/nu_pair_tophit_u_pair');

    &mcall(
           "RETRIEVE UNIQUE PAIRED READS",
           \&retrieve_status_from_psl,
           $temp . '/nu_pair_tophit_u_pair',
           $temp . '/nu_pair_tophit.psl',
           $temp . '/u_pair_from_nu_pair_tophit.psl'
          );

    &mcall(
           "FILTER INCONSISTENT PAIRED READS",
           \&filter_paired_reads,
           $temp . '/u_pair_from_nu_pair_tophit.psl',
           $temp . '/u_pair_from_nu_pair_tophit_consistent.psl',
           $temp . '/u_pair_from_nu_pair_tophit_inconsistent',
           $INTRON_MAX
          );

    my $comm_merge =
        'cat ' 
      . $temp
      . '/u_pair_consistent.psl '
      . $temp
      . '/u_pair_from_nu_pair_tophit_consistent.psl >'
      . $temp
      . '/u_pair_ALL.psl';
    mlog("== COPY ==", $comm_merge);
    system($comm_merge );

## find junctions

    my $junction_path = $current_path . '/junctions';
    mkdir($junction_path);
    &mcall(
           "FIND ALL JUNCTIONS",      \&find_junctions,
           $temp . '/u_pair_ALL.psl', $junction_path . '/junctions.junc',
           $INTRON_MIN,               $INTRON_MAX,
           $FLANKSIZE
          );

    &mcall("FIND UNIQUE JUNCTIONS",
           \&unique_junctions,
           $junction_path . '/junctions.junc',
           $junction_path . '/unique.junc');

    &mcall("FILTER JUNCTIONS",
           \&filter_human_splicesite,
           $junction_path . '/unique.junc',
           $junction_path . '/unique_GTAG.junc', $GENOME_DIR);

    my $end_time = time();
    mlog(" Finished. Total time: ", $end_time - $start_time, " seconds");

}    # END_OF_RUN
#################################### Other Subroutines  #################################

=head2   run_blat
 Title   : run_blat
 Usage   : 
 Function: run multiple Blat processes on a server
 Returns : 
 Args    : 
=cut

sub run_blat
{

    my ($fasta, $db, $chunk, $tilesize, $stepsize, $minscore, $rep, $cp) = @_;
    my $repmatch = int(1024.0 * $tilesize / $stepsize);
    my $parameters;

    if (not $rep)
    {
        $parameters =
          "-noHead  -minScore=$minscore -stepSize=$stepsize -tileSize=$tilesize";
    }
    else
    {
        $parameters =
          "-noHead -repMatch=$repmatch -minScore=$minscore -stepSize=$stepsize -tileSize=$tilesize";
    }

######  check the databases for BLAT
    unless (-f $db)
    {
        croak "Can't find BLAT database file as $db\n";
    }

## You must run SplicePL in the same folder as fasta file

    unless (-f $fasta)
    {
        croak "Can't find $fasta file in $cp. Exit!";
    }

    if (!-d $cp . '/scripts')
    {
        mkdir($cp . '/scripts');
    }

    for (my $i = 1 ; $i <= $chunk ; $i++)
    {
        my $command = <<EOF;
#!/bin/sh' 
date >> $cp/run.$i.log
echo blat $db $parameters $fasta.$i $fasta.$i.out >> $cp/run.$i.log
blat $db $parameters $fasta.$i $fasta.$i.out
date >> $cp/run.$i.log
EOF

        my $of = 'run.' . $i . '.sh';
        my $out_fh = FileHandle->new($cp . '/scripts/' . $of, "w");
        print $out_fh $command;
        undef $out_fh;
    }

    for my $p (1 .. $chunk)
    {
        my $pid = fork();
        if ($pid == -1)
        {
            die;
        }
        elsif ($pid == 0)
        {

            my $blat_process = "bash $cp/scripts/run.$p.sh ";
            exec $blat_process or die;
        }
    }
    while (wait() != -1) { }
    mlog("All $chunk BLAT jobs Done");

## merge BLAT output files into one file
    my $s = 'cat ';

    for my $j (1 .. $chunk)
    {
        $s .= $fasta . '.' . $j . '.out' . ' ';
    }

    my $psl_fn   = $cp . '/temp/pairedreads.psl';
    my $comm_cat = $s . '> ' . $psl_fn;
    mlog("== RUN BLAT ==", $comm_cat);
    system($comm_cat);

}

=head2 split_fasta
 Title   : split_fasta
 Usage   : 
 Function: split a fasta file into N pierces
 Returns : NULL. Side effect is number of chunks of files will be created
 Args    : Three scalars, first is the path to the fasta file,
           Second is the number of chunks
           Third is the current path
=cut

sub split_fasta
{
    my ($fasta, $chunk, $cp) = @_;
    my $fasta_no = 0;

    my $in = FileHandle->new($fasta, "r");
    while (<$in>)
    {
        $fasta_no++ if /^>/;
    }

    seek($in, 0, 0);

    my $chunk_size = int($fasta_no / $chunk);
    mlog("== SPLIT FASTA ==  No of fasta records ",
         $fasta_no, " to $chunk pieces each of $chunk_size");

    if ($chunk == 1)
    {
        system("cp " . $fasta . " " . $fasta . ".1");
        mlog("The fasta file is copied to $fasta" . '.1');
    }
    $chunk_size++ if $chunk_size % 2 == 1;

    foreach my $i (1 .. ($chunk - 1))
    {
        my $out = $fasta . "." . $i;
        my $out_fh = FileHandle->new($out, "w");
        foreach my $j (1 .. $chunk_size)
        {
            my $l = <$in>;
            print $out_fh $l;
            $l = <$in>;
            print $out_fh $l;
        }
        undef $out_fh;
    }

    my $out = $fasta . "." . $chunk;

    my $out_fh = FileHandle->new($out, "w");
    while (<$in>) { print $out_fh $_; }
    $in->close;
    $out_fh->close;
}

=head2 mark_reads_status
 Title   : mark_reads_status
 Usage   : mark_reads_status("/path/to/psl_file","/path/to/fasta_file")
 Function: Each read will be marked as
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

sub mark_reads_status
{
    my ($psl, $fasta) = @_;

    my $status              = $psl . '.status';
    my %pairedreads_status  = ();
    my %forwardreads_status = ();
    my %reversereads_status = ();

    my $psl_fh   = FileHandle->new($psl,    "r");
    my $out_fh   = FileHandle->new($status, "w");
    my $fasta_fh = FileHandle->new($fasta,  "r");

    my $counter       = 0;
    my $maxFastaID    = 0;    # in fasta file
    my $fastaTotalSeq = 0;    # total number of sequences in fasta file
    my $maxReadID     = 0;    # in psl file

## retrieve read sequence ID
## we assume the read ID is like 12345a, 12345b and it is CONTINUALLY numbered from 1 to maxReadID
    while (<$fasta_fh>)
    {
        if (/^>\D*(\d+)/)
        {
            $fastaTotalSeq++;
            if ($maxFastaID < $1)
            {
                $maxFastaID = $1;
            }
        }
    }

    mlog("== MARK READS STATUS ==",
         "total fasta sequences: ",
         $fastaTotalSeq, " maxFastaID: ", $maxFastaID);

    while (<$psl_fh>)
    {

        # col 10 is the read ID
        my $id = (split /\t/)[9];

        if ($id =~ /(\d+)(.)/)
        {

            #print $id, "\t", $1, "\t", $2, "\n";
            if ($1 > $maxReadID)
            {
                $maxReadID = $1;
            }
            if ($2 eq 'f')
            {
                $forwardreads_status{$1}++;
            }
            elsif ($2 eq 'r')
            {
                $reversereads_status{$1}++;
            }
            $pairedreads_status{$1}++;
        }
    }

    #print "maxReadID: $maxReadID \n";

    if ($maxReadID > $maxFastaID)
    {
        croak "maxReadID larger than maxFastaID. Exit!\n";
    }

    foreach my $key (1 .. $maxFastaID)
    {
        if (!$pairedreads_status{$key})
        {
            $pairedreads_status{$key} = "nomap";
        }
        print $out_fh $key, "\t", $pairedreads_status{$key}, "\t";

        if ($pairedreads_status{$key} == 1)
        {
            if ($forwardreads_status{$key} == 1)
            {
                $pairedreads_status{$key} = "u_forward";
            }
            else
            {
                $pairedreads_status{$key} = "u_reverse";
            }

        }
        else
        {

            if (   $forwardreads_status{$key} == 1
                && $reversereads_status{$key} == 1)
            {
                $pairedreads_status{$key} = "u_pair";
            }
            elsif (   $forwardreads_status{$key} >= 1
                   && $reversereads_status{$key} >= 1)
            {
                $pairedreads_status{$key} = "nu_pair";
            }
            elsif (   $forwardreads_status{$key} > 1
                   && $reversereads_status{$key} < 1)
            {
                $pairedreads_status{$key} = "nu_forward";
            }
            elsif (   $forwardreads_status{$key} < 1
                   && $reversereads_status{$key} > 1)
            {
                $pairedreads_status{$key} = "nu_reverse";
            }

        }

        print $out_fh $pairedreads_status{$key}, "\n";

    }

    undef $out_fh;
    undef $fasta_fh;
    undef $psl_fh;
}

### Find status, which is one of
###  u_pair, u_forward, u_reverse, nu_pair, nu_forward, nu_reverse, nomap
###

sub find_status
{
    my ($status, $psl, $out) = @_;

    my $status_psl = $psl . '.status';

    my $in     = FileHandle->new($status_psl, "r");
    my $out_fh = FileHandle->new($out,        "w");

    while (<$in>)
    {
        my $a = (split /\t/)[2];
        if ($a =~ /^$status/)
        {
            print $out_fh $_;
        }

    }

    undef $in;
    undef $out_fh;
}

sub retrieve_status_from_psl
{
### given a list of u_pair read ID such as 123
### get the blat result in the psl file
### psl file has no HEADERS by -noHead when running blat

    my ($status, $psl, $out) = @_;
    my %ids = ();

    my $status_fh = FileHandle->new($status, "r");
    my $psl_fh    = FileHandle->new($psl,    "r");
    my $out_fh    = FileHandle->new($out,    "w");

    while (<$status_fh>)
    {
        chomp;
        my $id = (split /\t/)[0];

        my $idf = $id . 'f';
        my $idr = $id . 'r';
        $ids{$idf}++;
        $ids{$idr}++;
    }

    while (<$psl_fh>)
    {
        my $i = (split /\t/)[9];
        if ($ids{$i})
        {
            print $out_fh $_;
        }
    }
    undef $psl_fh;
    undef $status_fh;
    undef $out_fh;
}

## exclude reads
##   1. mapped to different chromosome
##   2. mapped to same strand
##   3. too far away

sub filter_paired_reads
{
    my ($psl, $out1, $out2, $maxdist) = @_;

    my $psl_fh  = FileHandle->new($psl,  "r");
    my $out1_fh = FileHandle->new($out1, "w");
    my $out2_fh = FileHandle->new($out2, "w");

    my $count = 1;

    while (my $line1 = <$psl_fh>)
    {
        my ($for_id, $for_str, $for_chr, $for_start) =
          (split /\t/, $line1)[9, 8, 13, 15];

        my $m1;
        my $m2;

        my $line2 = <$psl_fh>;
        my ($rev_id, $rev_str, $rev_chr, $rev_start) =
          (split /\t/, $line2)[9, 8, 13, 15];

        if ($for_id =~ /(\d+)f/)
        {
            $m1 = $1;
        }
        else
        {
            croak $for_id,
              "does not match the format of forward ID (e.g., 1234f) \n";
        }

        if ($rev_id =~ /(\d+)r/)
        {
            $m2 = $1;
        }
        else
        {
            croak $rev_id,
              "does not match the format of reverse ID (e.g., 1234r) \n";
        }

        if ($m1 ne $m2)
        {
            croak "forward read ID does not match reverse read ID. Exit.\n";
        }

        if ($for_str eq $rev_str || $for_chr ne $rev_chr)
        {
            ## strand and chromosome
            if ($for_str eq $rev_str && $for_chr ne $rev_chr)
            {
                print $out2_fh $count, " same str: ",
                  join("\t", $for_id, $for_str, $rev_id, $rev_str), "\t";
                print $out2_fh $count++, " diff chr: ",
                  join("\t", $for_id, $for_chr, $rev_id, $rev_chr), "\n";
            }
            elsif ($for_str eq $rev_str)
            {
                print $out2_fh $count++, " same str: ",
                  join("\t", $for_id, $for_str, $rev_id, $rev_str), "\n";
            }
            elsif ($for_chr ne $rev_chr)
            {
                print $out2_fh $count++, " diff chr: ",
                  join("\t", $for_id, $for_chr, $rev_id, $rev_chr), "\n";
            }

        }
        else
        {
            ## distance
            if (abs($for_start - $rev_start) >= $maxdist)
            {
                print $out2_fh $count++, " too far: ",
                  join("\t", $for_id, $for_chr, $rev_id, $rev_chr), "\n";
            }
            else
            {
                print $out1_fh $line1, $line2;
            }
        }
    }
    undef $out1_fh;
    undef $out2_fh;
    undef $psl_fh;
}

sub filter_nu_reads_from_psl
{
    my ($psl, $out1, $out2) = @_;
    my %h;
    my %count_max;
    my $previous;
    my $maxline;

    my $psl_fh  = FileHandle->new($psl,  "r");
    my $out1_fh = FileHandle->new($out1, "w");
    my $out2_fh = FileHandle->new($out2, "w");

    my @a = split /\t/, <$psl_fh>;
    $h{$a[9]} = $a[0];    ## col 10 is the seqID  ## col1 is the score
    $previous = $a[9];

    while (<$psl_fh>)
    {
        my @a = split /\t/;
        if (!$h{$a[9]})
        {
            if ($count_max{$previous})
            {
                print $out2_fh '?' . $previous, "\t", $count_max{$previous};
            }
            else
            {
                print $out1_fh $maxline;
            }
            $previous = $a[9];
        }
        if ($h{$a[9]} < $a[0])
        {
            $h{$a[9]}         = $a[0];
            $maxline          = $_;
            $count_max{$a[9]} = undef;
        }
        elsif ($h{$a[9]} == $a[0])
        {
            $count_max{$a[9]} = $maxline . $a[9] . "\t" . $_;
        }
    }
    undef $out1_fh;
    undef $out2_fh;
    undef $psl_fh;
}

=head2 find_junctions
Function : find junctions from a psl file 
 col  9: strand
 col 10: seqID
 col 11: qSize
 col 14: chr
 col 18: block no
 col 19: blockSize
 col 21: tStart  (target start position)
 col start from 1 (not zero)
=cut

sub find_junctions
{
    my ($psl, $out, $mingap, $maxgap, $minblocksize) = @_;

    my $first_row = 1;

    my $psl_fh = FileHandle->new($psl, "r");
    my $out_fh = FileHandle->new($out, "w");

    while (<$psl_fh>)
    {
        if ($first_row == 1)
        {
            my @a = split /\t/;
            if ($a[10] <= 60 && $a[10] >= 20)
            {
                $minblocksize = 10 if $minblocksize < 20;
            }
            elsif ($a[10] > 60)
            {    #query size >=60, might need larger blocksize
                $minblocksize = 20 if $minblocksize < 20;
            }
            else
            {
                croak "PSL file has abnormal qSize: $a[10].";
            }
            mlog("== FIND JUNCTIONS ==", "minblocksize=$minblocksize");
            $first_row = undef;
        }
        my @a = split /\t/;
        my @tstarts = split /,/, $a[20];

        my @blocksize = split /,/, $a[18];
        my $blockno = $a[17];
        if ($blockno > 1)
        {
            foreach (0 .. $#tstarts - 1)
            {
                my $gapsize =
                  $tstarts[$_ + 1] + 1 - ($tstarts[$_] + $blocksize[$_]);
                if (   $gapsize >= $mingap
                    && $gapsize <= $maxgap
                    && $blocksize[$_] >= $minblocksize
                    && $blocksize[$_ + 1] >= $minblocksize)
                {

                    print $out_fh $a[13], ":", $tstarts[$_] + $blocksize[$_],
                      "-",
                      $tstarts[$_ + 1] + 1, "\t", $a[9], "\t", $a[8], "\t",
                      $tstarts[$_ + 1] + 1 - ($tstarts[$_] + $blocksize[$_]),
                      "\n";
                }
            }
        }
    }
    undef $psl_fh;
    undef $out_fh;
}

=head2 unique_junctions 
Function : Find unique junctions and count the coverage for each junction
=cut

sub unique_junctions
{
    my ($in, $out) = @_;

    my %h;
    my %line;
    my %l;
    my $unique_junctions_no = 0;

    my $in_fh  = FileHandle->new($in,  "r");
    my $out_fh = FileHandle->new($out, "w");

    while (<$in_fh>)
    {
        if (/^chr/)
        {
            chomp;
            my @a = split /\t/;
            $h{$a[0]}++;
            $line{$a[0]} = $_;
            $l{$a[0]} .= '_' . $a[1];
        }
    }

    foreach (sort keys %h)
    {
        if ($h{$_})
        {
            print $out_fh $line{$_}, "\t", $l{$_}, "\n";
            $unique_junctions_no++;
        }
    }

    undef $in_fh;
    undef $out_fh;
    return $unique_junctions_no;
}

=head2 filter_human_splicesite
 Function:
    Filter the junctions with the splicing site patterns
    1. GT...AG
    2. GC...AG
    3. AT...AC
=cut

sub filter_human_splicesite
{
    my ($junc, $out, $dir) = @_;

    my ($splus, $sneg, %chr_seq);

    my $junc_fh = FileHandle->new($junc, "r");
    my $out_fh  = FileHandle->new($out,  "w");

    my @chr_names = (
                     1,  2,  3,   4,   5,  6,  7,  8,  9,  10,
                     11, 12, 13,  14,  15, 16, 17, 18, 19, 20,
                     21, 22, 'X', 'Y', 'M'
                    );

### Load each chromosome to memory, for human genome 19, it uses 3GB of RAM
    foreach my $chr (@chr_names)
    {
        my $file = $dir . '/chr' . $chr . '.fa';

        mlog("== FILTER SPLICE SITES == Loading chromosome $file");
        my $chr_fh = FileHandle->new($file, "r");
        my $id = <$chr_fh>;

        my $seq = '';
        while (<$chr_fh>)
        {
            chomp;
            $seq .= $_;
        }
        $chr_seq{'chr' . $chr} = $seq;
        undef $chr_fh;
    }

### We were not analyzing strand-specific RNA-seq data, so
### have to search both plus and minus strand for
### splicing site patterns

    my $good_junctions_no = 0;
    my $bad_junctions_no  = 0;
    while (<$junc_fh>)
    {
        chomp;
        my @a = split /\s+/;
        my ($chr, $start, $end) = split /[:-]/, $a[0];
        my $strand = $a[2];
        my $s;
        $splus =
            substr($chr_seq{$chr}, $start, 2)
          . substr($chr_seq{$chr}, $end - 3, 2);
        $sneg = &revcom($splus);

        if ($splus =~ /^GTAG$/i)
        {
            print $out_fh $_, "\tGOODSITEPLUS" . $strand . "GTAGSITE\n";
            $good_junctions_no++;
        }
        elsif ($sneg =~ /^GTAG$/i)
        {
            print $out_fh $_, "\tGOODSITENEG" . $strand . "GTAGSITE\n";
            $good_junctions_no++;
        }
        elsif ($splus =~ /^GCAG$/i)
        {
            print $out_fh $_, "\tGOODSITEPLUS" . $strand . "GCAGSITE\n";
            $good_junctions_no++;
        }
        elsif ($sneg =~ /^GCAG$/i)
        {
            print $out_fh $_, "\tGOODSITENEG" . $strand . "GCAGSITE\n";
            $good_junctions_no++;
        }
        elsif ($splus =~ /^ATAC$/i)
        {
            print $out_fh $_, "\tGOODSITEPLUS" . $strand . "ATACSITE\n";
            $good_junctions_no++;
        }
        elsif ($sneg =~ /^ATAC$/i)
        {
            print $out_fh $_, "\tGOODSITENEG" . $strand . "ATACSITE\n";
            $good_junctions_no++;
        }
        else
        {
            print $out_fh $_, "\tBADSITE$strand\n";
            $bad_junctions_no++;
        }
    }
    mlog(
        "== FIND HUMAN SPLICECITES == Good junctions number $good_junctions_no; Bad junctions number $bad_junctions_no"
    );
    undef $out_fh;
    undef $junc_fh;
    return $good_junctions_no;
}

=head2 revcom
 Function: reverse and complement a DNA sequence
=cut

sub revcom
{
    my ($seq) = @_;

    my $r = reverse($seq);
    $r =~ tr/atcgATCG/tagcTAGC/;

    return $r;

}

=head2 check_total_memory
 Function: Retrieves the memory installed on this machine
=cut

sub check_total_memory
{
    my ($os) = @_;

    sub memerror
    {
        mlog("Unable to determine total memory");
    }

    my ($physical_memory, $swap_memory, $duflags);
    $os = $os || `uname`;

    if ($os =~ /Linux/)
    {
        $physical_memory = `free -b | grep Mem | awk '{print \$2}'` or memerror;
    }
    elsif ($os =~ /Darwin/)
    {
        $physical_memory = `sysctl -n hw.memsize` or memerror;
    }
    elsif ($os =~ /NetBSD|OpenBSD/)
    {
        $physical_memory = `sysctl -n hw.physmem` or memerror;
        if ($physical_memory < 0)
        {
            $physical_memory = `sysctl -n hw.physmem64` or memerror;
        }

    }
    elsif ($os =~ /BSD/)
    {
        $physical_memory = `sysctl -n hw.realmem`;

    }
    elsif ($os =~ /SunOS/)
    {
        $physical_memory = `/usr/sbin/prtconf | grep Memory | cut -f 3 -d ' '`
          or memerror;
        chomp($physical_memory);
        $physical_memory = $physical_memory * 1024 * 1024;
    }
    else
    {
        $physical_memory = 0;
        mlog("Can't figure out the operating system.");

    }

    chomp($physical_memory);
    return $physical_memory;
}

=head2 usage
 Function: Print out usage message
=cut

sub usage
{
    print <<"END_USAGE";
Usage: $0 --forward=read1.fa --reverse=read2.fa --genome=/data/genome --genome_2bit=hg19.2bit --processes 2 

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

=head2 prepare_paired_reads
Title: prepare_paired_reads
Usage : prepare_paired_reads("/path/to/read1.fa", "/path/to/read2.fa", "/path/to/output.fa")
Function : Check fasta format; assign numerical ID such as 123f to forward, 123r to reverse read; combine two reads files to 1 file; make sure sequences are of same number in forward and reverse reads
Returns : 
Args : First string is the path to forward read file
        Second string is the path to reverse read file
        Third string is the path to output fasta file
=cut

sub prepare_paired_reads
{
    my ($forward_fn, $rev_fn, $out_fn) = @_;

    # check fasta sequence number in two files
    my $f = &fasta_dnaseq_no($forward_fn);
    my $r = &fasta_dnaseq_no($rev_fn);

    if ($f == 0)
    {
        mlog(" $forward_fn is not a valid fasta file.");
        return 0;
    }

    if ($r == 0)
    {
        mlog(" $rev_fn is not a valid fasta file.");
        return 0;
    }

    if ($f != $r)
    {
        mlog(
            " $forward_fn has $f sequences, $rev_fn has $r sequences, not equal."
        );
        return 0;
    }

    my $fh = FileHandle->new($forward_fn, "r");
    my $rh = FileHandle->new($rev_fn,     "r");
    my $oh = FileHandle->new($out_fn,     "w");

    # read one sequence at a time
    my $i     = 1;    #i sequence counter
    my $f_seq = '';
    my $r_seq = '';

    # discard first line
    <$fh>;
    <$rh>;

    # loop through forward read file
    while (my $f_nl = <$fh>)
    {
        while (defined $f_nl && $f_nl !~ /^>/ && $f_nl !~ /^$/)
        {
            chomp($f_nl);
            $f_seq .= $f_nl;
            $f_nl = <$fh>;
        }
        print $oh '>' . $i . "f\n" . $f_seq, "\n";
        $f_seq = '';

        # reverse read
        my $r_nl = <$rh>;
        while (defined $r_nl && $r_nl !~ /^>/ && $r_nl !~ /^$/)
        {
            chomp($r_nl);
            $r_seq .= $r_nl;
            $r_nl = <$rh>;
        }
        print $oh '>' . $i . "r\n" . $r_seq, "\n";
        $r_seq = '';
        $i++;
    }

    undef $fh;
    undef $rh;
    undef $oh;

    return $i - 1;
}

=head2 fasta_dnaseq_no
 Function: Find total number of sequences in a fasta file; Also check the sanity of fasta file
 Args: The path to the fasta file
=cut

sub fasta_dnaseq_no
{
    my ($filename) = @_;
    my $i = 0;
    my $fh = FileHandle->new($filename, "r");
    while (<$fh>)
    {
        if (/^>/)
        {
            $i++;

        }
        else
        {
            chomp;
            if (/^$/)
            {
                mlog("Find Empty line in $filename Line: $.");
                return 0;
            }

            if (/[^ATCGatcgNn]/)
            {
                mlog("Find Illegal character in $filename Line: $.");
                return 0;
            }
        }
    }
    undef $fh;
    return $i;
}

=head2 gettime
  Function: Return current time formatted as: Tue Jun 22 08:13:20 2010
=cut

sub gettime
{
    my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
    return $now_string;
}

=head2 mlog
 Function: Log time and message to stdout
=cut

sub mlog
{
    my $message = join " ", @_;
    print &gettime(), " ", $message, "\n";

}

=head2 mcall
 Function: Log subroutine name and arguments, then Call subroutine with arguments
=cut

sub mcall
{
    my $annotation = shift;
    my ($caller, @arg) = @_;
    &mlog(join(" ", "==", $annotation, "==", @arg));
    &{$caller}(@arg);
}

1;
