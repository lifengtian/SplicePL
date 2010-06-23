use Test::More tests => 2;

#1
BEGIN {
use_ok( 'SplicePL' );
}

is ( SplicePL::check_total_memory("Win32"), 0,'check_total_memory() Win32');
