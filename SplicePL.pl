#!/usr/bin/perl -w

### blat pipeline on pierce server
### Lifeng Tian
### Apr 2010

my $usage="perl blatpipeline.pl <fasta> <no. BLAT processes> <tileSize> <stepSize> <running folder> <minScore>\n";
if ( @ARGV != 6 ) {
	die $usage;
}

my $myhome="/data/share/LifengTian/projects/RNAseq";
my $pipeline=$myhome.'/pipeline';   # blat pipeline


my $fasta=$ARGV[0];  # input fasta file
my $chunk=$ARGV[1];  # number of threads to parallelize
my $tileSize=$ARGV[2];
my $stepSize=$ARGV[3];
my $cp=$ARGV[4];     # running folder where fasta file is in
my $minScore=$ARGV[5]; 

## blast parameters
my $repMatch = int(1024.0*$tileSize/$stepSize);
my $parameters="-noHead -repMatch=$repMatch -minScore=$minScore -stepSize=$stepSize -tileSize=$tileSize  -out=pslx";
#my $parameters="-noHead -repMatch=$repMatch -minScore=$minScore -stepSize=$stepSize -tileSize=$tileSize -ooc=$pipeline/hg19.${tileSize}ooc.masked -out=pslx";

unless ( -f $cp.'/'.$fasta) {
    die "Can't find $fasta file!";
}

mkdir( $cp.'/scripts');


### Generate shell scripts to run Blat
for ( my $i=1;$i<=$chunk;$i++) {
    my $comm=<<EOF;
#!/bin/sh' 
date >> $cp/run.$i.log
echo $pipeline/blat $pipeline/hg19.fa.masked.2bit $parameters $cp/$fasta.$i $cp/$fasta.$i.out >> $cp/run.$i.log
$pipeline/blat $pipeline/hg19.fa.masked.2bit $parameters $cp/$fasta.$i $cp/$fasta.$i.out
date >> $cp/run.$i.log
EOF

	my $of='run.'.$i.'.sh';
	open OUT, ">$cp/scripts/".$of or die "Error open $of for writing. ";
	print OUT $comm;
}

### Make pieces of input fasta file

breakup_fasta($fasta,$chunk,$cp);



### Confirm with user the pipeline status

### Run the pipeline
for my $p (1..$chunk) {
   my $pid = fork();
   if ($pid == -1) {
       die;
   } elsif ($pid == 0) {
    
    	my $comm="bash $cp/scripts/run.$p.sh ";
  	exec $comm or die;
   }
}
while (wait() != -1) {}
print "All $chunk jobs Done \n";




### Run analysis



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
    $outfilename = $cp.'/'.$fastafile. "." . $numpieces;
    open(OUTFILE, ">$outfilename");
    while(my $line = <INFILE>) {
        print OUTFILE $line;
    }
    close(OUTFILE);
    return 0;
}

