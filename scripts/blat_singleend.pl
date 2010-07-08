#!/usr/bin/perl 

# PERL script for running BLAT and find uniquely-mapped reads
#
# Please direct questions and support issues to <lifeng2@mail.med.upenn.edu>
#
# Copyright Lifeng Tian
#
# You may distribute this module under the same terms as PERL itself

# POD documentation - main docs before the code

=head1 NAME
	blat_singleend.pl Run BLAT and Report uniquely-mpped reads above minscore and minidentity
	
=head1 SYNOPSIS
	perl blat_singleend.pl  --fasta reads.fa --genome_dir=/data/share/db/human/hg19 --genome_2bit_filename=hg19.2bit --minScore=60 --processes=6 --minidentity=90
	
=head1 DESCRIPTION
       This script simply run BLAT on N processors by splitting the input fasta files to N pieces; then find
       uniquely mapped reads from BLAT output.

       Please select stringent threshold for minscore, which is the matches minus the 
               mismatches minus some sort of gap penalty. The rule of thumb is, e.g., if you want at least 90 bases
       matching the reference sequence, make --minscore=90.
       
       You can also assign --minidentity=90, requiring a 90% minimal identity across the whole read length.


=head1 PREREQUISITES  
	Before running the script, please make sure:
	1. Either BLAT or gfServer/gfClient binaries is installed. The most recent version (Jun 25, 2010) is 34x7 .
	2. A 2BIT format file is prepared from reference genome fasta file(s). 

=head1 TO DO
       nothing on the list yet.

=head1 FEEDBACK

=head1 AUTHOR - Lifeng Tian

Email: lifeng2@mail.med.upenn.edu
	
=cut 

# Let the code begin ...

use strict;
use Carp;
use Cwd 'abs_path';
use Env '@PATH';
use FileHandle;
use POSIX qw(strftime);

use Getopt::Long;

my $USAGE                = '';
my $VERSION              = '0.1';
my $GENOME_DIR           = undef;
my $GENOME_2BIT_FILENAME = undef;
my $PROCESSES            = undef;
my $TILESIZE             = 11;
my $STEPSIZE             = 5;
my $MINSCORE             = 30;
my $REPMATCH             = undef;
my $GFSERVER             = undef;
my $SERVER               = 'localhost';
my $PORT                 = 5000;
my $FASTA_INPUT          = undef;
my $MINIDENTITY          = 0;
   # MINIDENTITY*length_of_read must be smaller than or equal to the readScore. By default, MINIDENTITY is set as 0.

my $getopt = Getopt::Long::GetOptions(
	'help|usage'             => \$USAGE,
	'genome_dir=s'           => \$GENOME_DIR,
	'genome_2bit_filename=s' => \$GENOME_2BIT_FILENAME,
	'processes=i'            => \$PROCESSES,
	'tilesize=i'             => \$TILESIZE,
	'stepsize=i'             => \$STEPSIZE,
	'minscore=i'             => \$MINSCORE,
	'repmatch'               => \$REPMATCH,
	'gfServer'               => \$GFSERVER,
	'server=s'               => \$SERVER,
	'port=i'                 => \$PORT,
	'fasta=s'                => \$FASTA_INPUT,
	'minidentity=i'          => \$MINIDENTITY,
);

# commandline usage
&usage(), exit(1) if ($USAGE);

#
if (   not $GENOME_DIR
	or not $GENOME_2BIT_FILENAME
	or not $FASTA_INPUT )
{
	if ( not $GENOME_DIR ) {
		print "Set --genome_dir=/path/to/fasta .\n";
	}

	if ( not $GENOME_2BIT_FILENAME ) {
		print "Set --genome_2bit_filename=xxx.2bit \n";
	}

	if ( not $FASTA_INPUT ) {
		print "set --fasta=xxx.fa \n";
	}
	&usage();
	exit(1);
}

my $start_time = time();
&mlog("blat_singleend.pl($VERSION) started");

print "Parameters\n 
GENOME_DIR = $GENOME_DIR 
GENOME_2BIT_FILENAME = $GENOME_2BIT_FILENAME 
FASTA_INPUT = $FASTA_INPUT
PROCESSES = $PROCESSES
TILESIZE = $TILESIZE
STEPSIZE = $STEPSIZE
MINSCORE = $MINSCORE
REPMATCH = $REPMATCH 
GFSERVER = $GFSERVER
SERVER = $SERVER
PORT = $PORT
MINIDENTITY = $MINIDENTITY
\n";

my $error = undef;

### check genome_2bit size and total available physical memory

if ( -f $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME ) {
	my $filesize = -s $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME;
	mlog( "== MEMORY ==", $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME,
		" size ", int( $filesize / 1024 / 1024 ), "MB" );
	my $total_mem = int( &check_total_memory() / 1024 / 1024 / 1024 );
	my $one_process_mem;

	if ($GFSERVER) {
		$one_process_mem = $filesize * 2 / 1024 / 1024 / 1024;
	}
	else {
		$one_process_mem = $filesize * 6 / 1024 / 1024 / 1024;
	}

	mlog("== MEMORY == Total physical memory is $total_mem GB");
	mlog("== MEMORY == Each process needs $one_process_mem GB");

	if ( $total_mem <= $one_process_mem ) {
		mlog(
"== MEMORY == Total physical memory < $one_process_mem , not enough to load the genome. Exit"
		);
		exit(1);
	}
	if ( !$PROCESSES ) {
		$PROCESSES = int( $total_mem / $one_process_mem );
		mlog(
			"== MEMORY == $PROCESSES processes will be used in the analysis" );
	}
	mlog( "== MEMORY == Predicted memory usage is ",
		$PROCESSES * $one_process_mem, "GB" );

	if ( $total_mem <= $PROCESSES * $one_process_mem ) {
		mlog(
"== MEMORY == Total physical memory is not enough to run $PROCESSES processes simultaneously. 
    The maximum number of processes you can run is ",
			int( $total_mem / $one_process_mem )
		);
		exit(1);
	}

}
else {
	$error .=
	  "!!ERROR!! $GENOME_DIR.'/'.$GENOME_2BIT_FILENAME does NOT exist.\n";
}

if ($error) {
	mlog($error);
	exit(1);
}

### Create temporary folder for intermediate files

my $current_path    = abs_path();
my $temp            = $current_path . '/temp';
my $fasta_fn        = $current_path . '/' . $FASTA_INPUT;
my $psl_blat_fn     = $temp . '/reads_blat.psl';
my $psl_gfClient_fn = $temp . '/reads_gfClient.psl';
my $psl_fn          = undef;

if ( !-d $temp ) {
	mkdir($temp);
}
else {
	mlog("== CREATE TEMP FOLDER == $temp exists");
}

### Start the pipeline here

&mcall( "SPLIT FASTA", \&split_fasta, $fasta_fn, $PROCESSES, $current_path );

if ($GFSERVER) {
	### check executable for gfServer and gfClient
	my $gfServer = 'gfServer';
	my $gfClient = 'gfClient';

	unless ( grep -x "$_/$gfServer", @PATH ) {
		$error .= "Can't find gfServer executables in @PATH\n";
	}

	unless ( grep -x "$_/$gfClient", @PATH ) {
		$error .= "Can't find gfClient executables in @PATH\n";
	}

	$psl_fn = $psl_gfClient_fn;
	&mcall(
		"RUN_GFSERVER",                            \&run_gfServer,
		$GENOME_DIR . '/' . $GENOME_2BIT_FILENAME, $PROCESSES,
		$TILESIZE,                                 $STEPSIZE,
		$REPMATCH,                                 $SERVER,
		$PORT,                                     $current_path
	);
	&mcall(
		"RUN_GFCLIENT", \&run_gfClient,
		$fasta_fn,      $fasta_fn . ".gfclient.psl",
		$PROCESSES,     $MINSCORE,
		$SERVER,        $PORT,
		$current_path
	);
	&mcall(
		"MARK READS ALIGNMENT STATUS FROM GFCLIENT OUTPUT",
		\&mark_reads_status_for_singleend,
		$psl_fn, $fasta_fn
	);
}
else {
	### check executable for BLAT
	my $blat_bin = 'blat';
	unless ( grep -x "$_/$blat_bin", @PATH ) {
		$error .= "Can't find blat executables in @PATH\n";
	}

	$psl_fn = $psl_blat_fn;

	&mcall(
		"RUN BLAT", \&run_blat,
		$fasta_fn,  $GENOME_DIR . '/' . $GENOME_2BIT_FILENAME,
		$PROCESSES, $TILESIZE,
		$STEPSIZE,  $MINSCORE,
		$REPMATCH,  $current_path
	);
	&mcall(
		"MARK READS ALIGNMENT STATUS",
		\&mark_reads_status_for_singleend,
		$psl_fn, $fasta_fn
	);
}

## unique-reads
&mcall( "FIND UNIQUE PAIRED READS",
	\&find_status, 'unique', $psl_fn, $temp . '/unique' );
&mcall(
	"RETRIEVE UNIQUE ALIGNMENTS",
	\&retrieve_status_from_psl, $temp . '/unique',
	$psl_fn, $temp . '/unique.psl'
);

my $end_time = time();
mlog( " Finished. Total time: ", $end_time - $start_time, " seconds" );

# END_OF_RUN
#################################### Other Subroutines  #################################

=head2   run_blat
 Title   :  run_blat
 Usage   :  
 Function:  run multiple Blat processes on a server
 Returns :  
 Args    :  First string is the absolute path to the input fasta file
 			Second string is the absolute path to the input 2bit file

=cut

sub run_blat {

	my ( $fasta, $db_2bit, $chunk, $tilesize, $stepsize, $minscore, $rep, $cp )
	  = @_;
	my $repmatch = int( 1024.0 * $tilesize / $stepsize );
	my $parameters;

	if ( not $rep ) {
		$parameters =
"-noHead  -minScore=$minscore -stepSize=$stepsize -tileSize=$tilesize";
	}
	else {
		$parameters =
"-noHead -repMatch=$repmatch -minScore=$minscore -stepSize=$stepsize -tileSize=$tilesize";
	}

######  check the databases for BLAT
	unless ( -f $db_2bit ) {
		croak "Can't find BLAT database file as $db_2bit\n";
	}

## You must run SplicePL in the same folder as fasta file

	unless ( -f $fasta ) {
		croak "Can't find $fasta file. Exit!";
	}

	if ( !-d $cp . '/scripts' ) {
		mkdir( $cp . '/scripts' );
	}

	for ( my $i = 1 ; $i <= $chunk ; $i++ ) {
		my $command = <<TEMPLATE;
#!/bin/sh' 
date >> $cp/run.$i.log
echo blat $db_2bit $parameters $fasta.$i $fasta.$i.out >> $cp/run.$i.log
blat $db_2bit $parameters $fasta.$i $fasta.$i.out
date >> $cp/run.$i.log
TEMPLATE

		my $of = 'run.' . $i . '.sh';
		my $out_fh = FileHandle->new( $cp . '/scripts/' . $of, "w" );
		print $out_fh $command;
		undef $out_fh;
	}

# spawn subprocesses and wait till each of them finish
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
	mlog("All $chunk BLAT jobs Done");

## merge BLAT output files into one file
	my $tmp = 'cat ';

	for my $j ( 1 .. $chunk ) {
		$tmp .= $fasta . '.' . $j . '.out' . ' ';
	}

	my $psl_fn   = $cp . '/temp/reads_blat.psl';
	my $comm_cat = $tmp . '> ' . $psl_fn;
	mlog( "== RUN BLAT ==", $comm_cat );
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

sub split_fasta {
	my ( $fasta, $chunk, $cp ) = @_;
	my $fasta_no = 0;

	my $in = FileHandle->new( $fasta, "r" );
	while (<$in>) {
		$fasta_no++ if /^>/;
	}

	seek( $in, 0, 0 );

	my $chunk_size = int( $fasta_no / $chunk );
	mlog( "== SPLIT FASTA ==  No of fasta records ",
		$fasta_no, " to $chunk pieces each of $chunk_size" );

	if ( $chunk == 1 ) {
		system( "cp " . $fasta . " " . $fasta . ".1" );
		mlog( "The fasta file is copied to $fasta" . '.1' );
	}
	$chunk_size++ if $chunk_size % 2 == 1;

	foreach my $i ( 1 .. ( $chunk - 1 ) ) {
		my $out = $fasta . "." . $i;
		my $out_fh = FileHandle->new( $out, "w" );
		foreach my $j ( 1 .. $chunk_size ) {
			my $l = <$in>;
			print $out_fh $l;
			$l = <$in>;
			print $out_fh $l;
		}
		undef $out_fh;
	}

	my $out = $fasta . "." . $chunk;

	my $out_fh = FileHandle->new( $out, "w" );
	while (<$in>) { print $out_fh $_; }
	$in->close;
	$out_fh->close;
}

=head2 mark_reads_status_for_singleend
 Title   : mark_reads_status_for_singleend
 Usage   : mark_reads_status_for_singleend("/path/to/psl_file","/path/to/fasta_file")
 Function: read will be marked as
            unique
            nonuniq
 Returns : 
=cut

sub mark_reads_status_for_singleend {
	my ( $psl, $fasta ) = @_;

	my $status       = $psl . '.status';
	my %reads_status = ();

	my $psl_fh   = FileHandle->new( $psl,    "r" );
	my $out_fh   = FileHandle->new( $status, "w" );
	my $fasta_fh = FileHandle->new( $fasta,  "r" );

	my $fastaTotalSeq = 0;    # total number of sequences in fasta file

	while (<$fasta_fh>) {
		$fastaTotalSeq++;
	}

	mlog( "== MARK READS STATUS ==", "total fasta sequences: ",
		$fastaTotalSeq );

	while (<$psl_fh>) {
		my $readLength = ( split /\t/ )[10];
		my $readScore  = ( split /\t/ )[0];
		if ( $readScore >= $readLength * $MINIDENTITY / 100 ) {

			# col 10 is the read ID
			my $id = ( split /\t/ )[9];
			$reads_status{$id}++;
		}
	}

	foreach my $key ( keys %reads_status ) {
		print $out_fh $key, "\t", $reads_status{$key}, "\t";

		if ( $reads_status{$key} == 1 ) {
			$reads_status{$key} = "unique";
		}

		else {
			$reads_status{$key} = "nonuniq";
		}

		print $out_fh $reads_status{$key}, "\n";

	}

	undef $out_fh;
	undef $fasta_fh;
	undef $psl_fh;
}

### Find status, which is one of
###  unique, nonuniq

sub find_status {
	my ( $status, $psl, $out ) = @_;

	my $status_psl = $psl . '.status';

	my $in     = FileHandle->new( $status_psl, "r" );
	my $out_fh = FileHandle->new( $out,        "w" );

	while (<$in>) {
		my $a = ( split /\t/ )[2];
		if ( $a =~ /^$status/ ) {
			print $out_fh $_;
		}

	}

	undef $in;
	undef $out_fh;
}

sub retrieve_status_from_psl {
### given a list of  read ID such as 123
### retrieve the blat result from the psl file
### psl file has no HEADERS by -noHead when running blat

	my ( $status, $psl, $out ) = @_;
	my %ids = ();

	my $status_fh = FileHandle->new( $status, "r" );
	my $psl_fh    = FileHandle->new( $psl,    "r" );
	my $out_fh    = FileHandle->new( $out,    "w" );

	while (<$status_fh>) {
		chomp;
		my $id = ( split /\t/ )[0];
		$ids{$id}++;
	}

	while (<$psl_fh>)
	{

		my $readLength = ( split /\t/ )[10];
		my $readScore  = ( split /\t/ )[0];
		if ( $readScore >= $readLength * $MINIDENTITY / 100 ) {

			# col 10 is the read ID
			my $id = ( split /\t/ )[9];
			if ( $ids{$id} ) {
				print $out_fh $_;
			}
		}
	}
	undef $psl_fh;
	undef $status_fh;
	undef $out_fh;
}

=head2 check_total_memory
 Function: Retrieves the memory installed on this machine
=cut

sub check_total_memory {
	my ($os) = @_;

	sub memerror {
		mlog("Unable to determine total memory");
	}

	my ( $physical_memory, $swap_memory, $duflags );
	$os = $os || `uname`;

	if ( $os =~ /Linux/ ) {
		$physical_memory = `free -b | grep Mem | awk '{print \$2}'` or memerror;
	}
	elsif ( $os =~ /Darwin/ ) {
		$physical_memory = `sysctl -n hw.memsize` or memerror;
	}
	elsif ( $os =~ /NetBSD|OpenBSD/ ) {
		$physical_memory = `sysctl -n hw.physmem` or memerror;
		if ( $physical_memory < 0 ) {
			$physical_memory = `sysctl -n hw.physmem64` or memerror;
		}

	}
	elsif ( $os =~ /BSD/ ) {
		$physical_memory = `sysctl -n hw.realmem`;

	}
	elsif ( $os =~ /SunOS/ ) {
		$physical_memory = `/usr/sbin/prtconf | grep Memory | cut -f 3 -d ' '`
		  or memerror;
		chomp($physical_memory);
		$physical_memory = $physical_memory * 1024 * 1024;
	}
	else {
		$physical_memory = 0;
		mlog("Can't figure out the operating system.");

	}

	chomp($physical_memory);
	return $physical_memory;
}

=head2 usage
 Function: Print out usage message
=cut

sub usage {
	print <<"END_USAGE";
Usage: $0 --fasta=read.fa --genome_dir=/data/genome --genome_2bit_filename=hg19.2bit --processes 2 --minidentity=90

--help                  Shows this help message
--fasta=name          Name of read file, fasta file only
--genome_dir=dir            Directory of genome sequences 
--genome_2bit_filename=name      Name of genome sequence file in 2bit format 
--processes             Number of processes for BLAT (default is 1)
--tilesize              Number of tileSize (default is 11)
--stepsize              Number of stepSize (default is 5)
--minscore              Number of minScore (default is 30)
--minidentity           Number of minIdentity (default is 0) WARNING! Use minidentity if you require minscore >= readlength * minidentity / 100
END_USAGE
}

=head2 fasta_dnaseq_no
 Function: Find total number of sequences in a fasta file; Also check the sanity of fasta file
 Args: The path to the fasta file
=cut

sub fasta_dnaseq_no {
	my ($filename) = @_;
	my $i = 0;
	my $fh = FileHandle->new( $filename, "r" );
	while (<$fh>) {
		if (/^>/) {
			$i++;

		}
		else {
			chomp;
			if (/^$/) {
				mlog("Find Empty line in $filename Line: $.");
				return 0;
			}

			if (/[^ATCGatcgNn]/) {
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

sub gettime {
	my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
	return $now_string;
}

=head2 mlog
 Function: Log time and message to stdout
=cut

sub mlog {
	my $message = join " ", @_;
	print &gettime(), " ", $message, "\n";

}

=head2 mcall
 Function: Log subroutine name and arguments, then Call subroutine with arguments
=cut

sub mcall {
	my $annotation = shift;
	my ( $caller, @arg ) = @_;
	&mlog( join( " ", "==", $annotation, "==", @arg ) );
	&{$caller}(@arg);
}

=head2   run_gfServer
 Title   : run_gfServer
 Usage   : 
 Function: run multiple gfServer processes on a server
 Returns : 
 Args    : 
=cut

sub run_gfServer {

	my ( $db_2bit, $chunk, $tilesize, $stepsize, $rep, $server, $port, $cp ) =
	  @_;
	my $repmatch = int( 1024.0 * $tilesize / $stepsize );
	my $parameters;

	if ( not $rep ) {
		$parameters = "-stepSize=$stepsize -tileSize=$tilesize -log";

	}
	else {
		$parameters =
		  "-repMatch=$repmatch  -stepSize=$stepsize -tileSize=$tilesize -log";
	}

######  check the databases for BLAT

	if ( !-d $cp . '/scripts' ) {
		mkdir( $cp . '/scripts' );
	}

	$port   = $port   || 5000;
	$server = $server || 'localhost';

	$db_2bit =~ /(.+)\.2bit/;
	my $db = $1 if defined $1;

	for ( my $i = 1 ; $i <= $chunk ; $i++ ) {
		unless ( -f $db . '.' . $i . '.2bit' ) {
			croak "Can't find BLAT database files\n";
		}

		my $command = <<TEMPLATE;
#!/bin/sh' 
date >> $cp/setup_gfServer.$i.log
echo gfServer start $server $port $db.$i.2bit $parameters=gfServer.log.$i >> $cp/setup_gfServer.$i.log
gfServer start $server $port $db.$i.2bit $parameters=gfServer.log.$i &
date >> $cp/setup_gfServer.$i.log
TEMPLATE

		my $of = 'setup_gfServer.' . $i . '.sh';
		my $out_fh = FileHandle->new( $cp . '/scripts/' . $of, "w" );
		print $out_fh $command;
		undef $out_fh;
		$port++;
	}

	for my $p ( 1 .. $chunk ) {
		my $pid = fork();
		if ( $pid == -1 ) {
			die;
		}
		elsif ( $pid == 0 ) {

			my $gfServer_process = "bash $cp/scripts/setup_gfServer.$p.sh ";
			exec $gfServer_process or die;
		}
	}
	while ( wait() != -1 ) { }
	mlog("All $chunk gfServer are up and running.");
	sleep(300);

	return 1;
}

=head2   run_gfClient
 Title   : run_gfClient
 Usage   : 
 Function: run multiple gfClient processes on a server
 Returns : 
 Args    : 
=cut

sub run_gfClient {

	my ( $input_fasta, $output_psl, $chunk, $minscore, $server, $port, $cp ) =
	  @_;

	my $parameters;

	$parameters = "-minScore=$minscore";

	if ( !-d $cp . '/scripts' ) {
		croak "scripts folder does not exist. Exit\n";
	}

	$port   = $port   || 5000;
	$server = $server || 'localhost';

	for ( my $i = 1 ; $i <= $chunk ; $i++ ) {

		my $command = <<TEMPLATE;
#!/bin/sh' 
date >> $cp/setup_gfClient.$i.log
echo  gfClient  $server $port / $input_fasta.$i $output_psl.$i $parameters >> $cp/setup_gfClient.$i.log
gfClient  $server $port / $input_fasta.$i $output_psl.$i $parameters
date >> $cp/setup_gfClient.$i.log
TEMPLATE

		my $of = 'setup_gfClient.' . $i . '.sh';
		my $out_fh = FileHandle->new( $cp . '/scripts/' . $of, "w" );
		print $out_fh $command;
		undef $out_fh;
		$port++;
	}

	for my $p ( 1 .. $chunk ) {
		my $pid = fork();
		if ( $pid == -1 ) {
			die;
		}
		elsif ( $pid == 0 ) {

			my $gfClient_process = "bash $cp/scripts/setup_gfClient.$p.sh ";
			exec $gfClient_process or die;
		}
	}
	while ( wait() != -1 ) { }
	mlog("All $chunk gfClient are Done.");
	system("killall gfServer");

### merge BLAT output files into one file
	my $tmp = 'cat ';

	for my $j ( 1 .. $chunk ) {
		$tmp .= $output_psl . '.' . $j . ' ';
	}

	my $psl_fn   = $cp . '/temp/reads_gfClient.psl';
	my $comm_cat = $tmp . '> ' . $psl_fn;
	mlog( "== RUN BLAT ==", $comm_cat );
	system($comm_cat);

}

1;
