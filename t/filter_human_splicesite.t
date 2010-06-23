use Test::More tests => 2;

BEGIN {
use_ok( 'SplicePL' );

}

my @chr_names = (
        1,  2,  3,  4,  5,  6,  7,  8,  9,   10,  11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y', 'M'
    );
foreach my $chr (@chr_names){
        print("chr$chr does not exist in $i.\n") unless -f '/data/share/db/human/hg19.fa.masked/chr'.$chr.'.fa' ;
    }       
ok ( SplicePL::filter_human_splicesite("t/junction12","t/junction12.out","/data/share/db/human/hg19.fa.masked")==4,'filter_human_splicesite() ');
