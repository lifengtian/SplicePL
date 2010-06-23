use Test::More tests => 5;

BEGIN { 
    use lib '..';

use_ok( 'SplicePL' );
}

is ( SplicePL::fasta_dnaseq_no("t/test3.fa"), 2,   'fasta_dnaseq_no() simple');
is ( SplicePL::fasta_dnaseq_no("t/test3_1.fa"), 0, 'fasta_dnaseq_no() non ATCG');
is ( SplicePL::fasta_dnaseq_no("t/test3_2.fa"), 2, 'fasta_dnaseq_no() lowercase, multiple lines');
is ( SplicePL::fasta_dnaseq_no("t/test3_3.fa"), 0, 'fasta_dnaseq_no() empty line');
