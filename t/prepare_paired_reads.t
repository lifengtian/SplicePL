use Test::More tests => 5;

BEGIN { 
	use_ok('SplicePL');
}


is ( SplicePL::prepare_paired_reads("t/test7_f.fa","t/test7_r.fa","t/test7_out.fa"),4,'prepare_paired_reads() check output line number');
is ( SplicePL::prepare_paired_reads("t/test3_1.fa","t/test7_r.fa","t/test7_out.fa"),0,'prepare_paired_reads() non ATCG');
is ( SplicePL::prepare_paired_reads("t/test3_2.fa","t/test7_r.fa","t/test7_out.fa"),0,'prepare_paired_reads() lowercase, multiple lines');
is ( SplicePL::prepare_paired_reads("t/test3_3.fa","t/test7_r.fa","t/test7_out.fa"),0,'prepare_paired_reads() empty line');
