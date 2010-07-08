use Test::More tests => 1;

`scripts/junction_to_bed.pl t/junction12 t/junction12.bed.new`;
my $out = `diff t/junction12.bed.new t/junction12.bed`;
chomp $out;

is $out, '' ,   'junction_to_bed';

