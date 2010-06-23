use Test::More tests => 3;

BEGIN { 
    use lib '..';
	
}

use_ok ('SplicePL');

is SplicePL::revcom("AAA"),'TTT','revcom';
is SplicePL::fasta_dnaseq_no("t/test3.fa"), 2,   'fasta_dnaseq_no() simple';

